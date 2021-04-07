import pytest
import networkx as nx

import os
import subprocess
import csv
import math
from copy import deepcopy
from glob import glob
import sys

from titan.parse_params import ObjMap, create_params
from titan.model import TITAN

# overwrite
@pytest.fixture
def params_integration(tmpdir):
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "simple_integration.yml"
    )
    return create_params(None, param_file, tmpdir)


@pytest.fixture
def make_model_integration(params_integration):
    def _make_model_integration():
        return TITAN(params_integration)

    return _make_model_integration


@pytest.mark.integration_deterministic
def test_model_runs():
    f = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "titan", "run_titan.py"
    )
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )

    subprocess.check_call([sys.executable, f, f"-p {param_file}"])
    assert True


@pytest.mark.integration_deterministic
def test_model_runs_sweep():
    f = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "titan", "run_titan.py"
    )
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )

    subprocess.check_call(
        [sys.executable, f, f"-p {param_file}", "-w model.seed.run:1:3"]
    )
    assert True


@pytest.mark.integration_deterministic
def test_model_reproducible(tmpdir):
    path_a = tmpdir.mkdir("result_a")
    path_b = tmpdir.mkdir("result_b")
    f = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "titan", "run_titan.py"
    )
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic_seeded.yml"
    )

    subprocess.check_call([sys.executable, f, f"-p {param_file}", f"-o {path_a}"])
    subprocess.check_call([sys.executable, f, f"-p {param_file}", f"-o {path_b}"])

    result_file_a = os.path.join(path_a, "basicReport.txt")
    result_file_b = os.path.join(path_b, "basicReport.txt")
    assert os.path.isfile(result_file_a)
    with open(result_file_a, newline="") as fa, open(result_file_b, newline="") as fb:
        reader_a = csv.DictReader(fa, delimiter="\t")
        res_a = []
        for row in reader_a:
            res_a.append(row)

        reader_b = csv.DictReader(fb, delimiter="\t")
        res_b = []
        for row in reader_b:
            res_b.append(row)

    for i in range(len(res_a)):
        assert res_a[i]["t"] == res_b[i]["t"]
        assert res_a[i]["run_id"] != res_b[i]["run_id"]
        assert res_a[i]["rseed"] == res_b[i]["rseed"]
        assert res_a[i]["pseed"] == res_b[i]["pseed"]
        assert res_a[i]["agents"] == res_b[i]["agents"]
        assert res_a[i]["hiv"] == res_b[i]["hiv"]
        assert res_a[i]["prep"] == res_b[i]["prep"]
        assert res_a[i]["deaths"] == res_b[i]["deaths"]
        assert res_a[i]["hiv_aids"] == res_b[i]["hiv_aids"]


@pytest.mark.integration_deterministic
def test_model_pop_write_read(tmpdir):
    path_a = tmpdir.mkdir("a")
    f = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "titan", "run_titan.py"
    )
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )

    subprocess.check_call(
        [sys.executable, f, "-p", param_file, "-o", path_a, "--savepop"]
    )

    saved_pop_path = glob(os.path.join(path_a, "pop", "*_pop.tar.gz"))[0]

    path_b = tmpdir.mkdir("b")
    subprocess.check_call(
        [sys.executable, f, "-p", param_file, "-o", path_b, "--poppath", saved_pop_path]
    )

    assert True


@pytest.mark.integration_deterministic
def test_model_settings_run(tmpdir):
    f = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "titan", "run_titan.py"
    )
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "integration_base.yml"
    )

    for item in os.listdir(
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "titan", "settings"
        )
    ):
        if "__" not in item and item != "base":
            path = tmpdir.mkdir(item)
            print(f"-----------Starting run for {item}-----------")
            subprocess.check_call(
                [
                    sys.executable,
                    f,
                    f"-S {item}",
                    f"-p {param_file}",
                    f"-o {path}",
                    "-e",
                ]
            )
            assert True


@pytest.mark.integration_deterministic
def test_agent_pop_stable_setting(tmpdir):
    """
    Agent population count stable across time
    """
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "integration_base.yml"
    )

    for item in os.listdir(
        os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "..", "titan", "settings"
        )
    ):
        if "__" not in item and item != "base":
            path = tmpdir.mkdir(item)
            os.mkdir(os.path.join(path, "network"))
            print(f"-----------Starting run for {item}-----------")
            params = create_params(item, param_file, tmpdir)
            model = TITAN(params)

            orig_pop = model.pop.all_agents.num_members()
            assert orig_pop == params.model.num_pop

            model.run(path)

            assert model.pop.all_agents.num_members() == orig_pop


@pytest.mark.integration_stochastic
def test_target_partners(make_model_integration, tmpdir):
    """
    If we increase the number of target partners, does the number of actual partners increase?
    """
    model_a = make_model_integration()
    model_a.params.outputs.network.edge_list = True

    path_a = tmpdir.mkdir("a")
    path_a.mkdir("network")
    path_b = tmpdir.mkdir("b")
    path_b.mkdir("network")

    # run with default bins (0-9)
    run_id_a = model_a.id
    model_a.run(path_a)

    # change the partner distribution mean upward for creating model b
    for bond in model_a.params.classes.bond_types:
        for race in model_a.params.classes.races:
            model_a.params.demographics[race].sex_type.MSM.drug_type[
                "None"
            ].num_partners[bond].vars[1].value *= 10
            model_a.params.demographics[race].sex_type.MSM.drug_type[
                "Inj"
            ].num_partners[bond].vars[1].value *= 10
    model_a.params.model.seed.run = model_a.run_seed
    model_a.params.model.seed.ppl = model_a.pop.pop_seed

    model_b = TITAN(model_a.params)
    run_id_b = model_b.id

    model_b.run(path_b)

    g_a_0 = nx.read_edgelist(
        os.path.join(path_a, "network", f"{run_id_a}_Edgelist_t0.txt"),
        delimiter="|",
        data=False,
    )
    g_a_10 = nx.read_edgelist(
        os.path.join(path_a, "network", f"{run_id_a}_Edgelist_t10.txt"),
        delimiter="|",
        data=False,
    )
    g_b_0 = nx.read_edgelist(
        os.path.join(path_b, "network", f"{run_id_b}_Edgelist_t0.txt"),
        delimiter="|",
        data=False,
    )
    g_b_10 = nx.read_edgelist(
        os.path.join(path_b, "network", f"{run_id_b}_Edgelist_t10.txt"),
        delimiter="|",
        data=False,
    )

    # should be at least 2x bigger
    assert (g_a_0.number_of_edges() * 1.5) < g_b_0.number_of_edges()
    assert (g_a_10.number_of_edges() * 1.5) < g_b_10.number_of_edges()


@pytest.mark.integration_stochastic
def test_prep_coverage(make_model_integration, tmpdir):
    """
    If we increase the target of prep coverage, does the incidence of hiv decrease?
    """
    path_a = tmpdir.mkdir("a")
    path_a.mkdir("network")
    path_b = tmpdir.mkdir("b")
    path_b.mkdir("network")

    model_a = make_model_integration()
    model_a.run(path_a)

    # change the coverage upward for creating model b, use same seeds
    model_a.params.prep.cap = 0.9
    model_a.params.prep.init = 0.9
    model_a.params.model.seed.run = model_a.run_seed
    model_a.params.model.seed.ppl = model_a.pop.pop_seed

    model_b = TITAN(model_a.params)
    model_b.run(path_b)
    print(
        f"model b prep world cap: {model_b.pop.geography.locations['world'].params.prep.cap}"
    )

    result_file_a = os.path.join(path_a, "basicReport.txt")
    result_file_b = os.path.join(path_b, "basicReport.txt")
    assert os.path.isfile(result_file_a)
    with open(result_file_a, newline="") as fa, open(result_file_b, newline="") as fb:
        reader_a = csv.DictReader(fa, delimiter="\t")
        res_a = {}
        for row in reader_a:
            if row["t"] not in res_a:
                res_a[row["t"]] = {"hiv": 0, "prep": 0}
            res_a[row["t"]]["hiv"] += float(row["hiv"])
            res_a[row["t"]]["prep"] += float(row["hiv"])

        reader_b = csv.DictReader(fb, delimiter="\t")
        res_b = {}
        for row in reader_b:
            if row["t"] not in res_b:
                res_b[row["t"]] = {"hiv": 0, "prep": 0}
            res_b[row["t"]]["hiv"] += float(row["hiv"])
            res_b[row["t"]]["prep"] += float(row["prep"])

    # at start, should be similar
    t0_hiv_a = res_a["0"]["hiv"]
    t0_hiv_b = res_b["0"]["hiv"]
    t0_diff = t0_hiv_a - t0_hiv_b
    assert math.isclose(t0_hiv_a, t0_hiv_b, abs_tol=15)  # within 15 agents
    assert res_a["0"]["prep"] < res_b["0"]["prep"]

    # at end, should be different
    t10_hiv_a = res_a["10"]["hiv"]
    t10_hiv_b = res_b["10"]["hiv"]
    t10_diff = t10_hiv_a - t10_hiv_b  # a should be higher
    assert res_a["10"]["prep"] < res_b["10"]["prep"]
    assert t10_hiv_a > t0_hiv_a
    assert t10_diff > t0_diff


@pytest.mark.integration_stochastic
def test_syringe_services(params_integration, tmpdir):
    """
    If we use syringe services, does the incidence of hiv decrease?
    """
    for race in params_integration.classes.races:
        params_integration.demographics[race].sex_type.MSM.drug_type.Inj.ppl = 1.0
        params_integration.demographics[race].sex_type.MSM.drug_type["None"].ppl = 0.0

    params_integration.model.num_pop = 500
    model_a = TITAN(params_integration)
    model_a.params.partnership.sex.frequency.Sex = (
        ObjMap(  # turn off sex to show only injection effects
            {"type": "bins", "bins": {1: {"prob": 1.0, "min": 0, "max": 1}}}
        )
    )
    path_a = tmpdir.mkdir("a")
    path_a.mkdir("network")
    path_b = tmpdir.mkdir("b")
    path_b.mkdir("network")

    # run with default bins (0-9)
    model_a.run(path_a)
    num_bonds = 0
    num_ag = 0
    for ag in model_a.pop.all_agents:
        if ag.hiv.active:
            num_ag += 1
            for ptnr in ag.get_partners(["Inj", "SexInj"]):
                if not ptnr.hiv.active:
                    num_bonds += 1
    assert num_bonds  # make sure there are serodiscordant partnerships

    # change the coverage upward for creating model b, use same seeds
    model_a.params.features.syringe_services = True
    model_a.params.model.seed.run = model_a.run_seed
    model_a.params.model.seed.ppl = model_a.pop.pop_seed

    model_b = TITAN(model_a.params)

    model_b.run(path_b)

    result_file_a = os.path.join(path_a, "basicReport.txt")
    result_file_b = os.path.join(path_b, "basicReport.txt")
    assert os.path.isfile(result_file_a)
    with open(result_file_a, newline="") as fa, open(result_file_b, newline="") as fb:
        reader_a = csv.DictReader(fa, delimiter="\t")
        res_a = {}
        for row in reader_a:
            if row["drug_type"] == "Inj":
                res_a[row["t"]] = res_a.get(row["t"], 0) + float(row["hiv"])

        reader_b = csv.DictReader(fb, delimiter="\t")
        res_b = {}
        for row in reader_b:
            if row["drug_type"] == "Inj":
                res_b[row["t"]] = res_b.get(row["t"], 0) + float(row["hiv"])

    # at start, should be similar
    t0_hiv_a = res_a["0"]
    t0_hiv_b = res_b["0"]
    t0_diff = t0_hiv_a - t0_hiv_b
    assert math.isclose(t0_hiv_a, t0_hiv_b, abs_tol=15)  # within 15 agents

    # at end, should be different
    t10_hiv_a = res_a["10"]
    t10_hiv_b = res_b["10"]
    t10_diff = t10_hiv_a - t10_hiv_b  # a should be higher
    assert t10_hiv_a > t0_hiv_a
    assert t10_hiv_b > t0_hiv_b
    assert t10_diff > t0_diff


@pytest.mark.integration_deterministic
def test_static_network(make_model_integration, tmpdir):
    model = make_model_integration()
    orig_rel_ids = [rel.id for rel in model.pop.relationships]

    # turn on static network - only affects run time, so fine have false during init
    model.params.features.static_network = True

    tmpdir.mkdir("network")

    for t in range(1, 10):
        model.time = t
        model.step(tmpdir)
        model.reset_trackers()

        curr_rel_ids = [rel.id for rel in model.pop.relationships]

        for rel in curr_rel_ids:
            assert rel in orig_rel_ids

        for rel in orig_rel_ids:
            assert rel in curr_rel_ids


@pytest.mark.integration_deterministic
def test_incar(params_integration, tmpdir):
    # turn on incar - initi is set to 0, so for these purposes, just run time
    params_integration.features.incar = True
    model = TITAN(params_integration)

    tmpdir.mkdir("network")

    # make sure we're starting from a clean slate
    init_incar = sum([1 for agent in model.pop.all_agents if agent.incar.active])
    assert init_incar == 0

    model.time = 1
    model.step(tmpdir)
    model.reset_trackers()

    time_1_incars = [agent for agent in model.pop.all_agents if agent.incar.active]
    time_1_incar_hiv_neg = [agent for agent in time_1_incars if not agent.hiv.active]
    should_release_t2 = [
        agent for agent in time_1_incars if agent.incar.release_time == 2
    ]

    assert len(time_1_incars) > 0
    assert len(time_1_incar_hiv_neg) > 0
    assert len(should_release_t2) > 0

    for agent in time_1_incars:
        assert agent.incar.time == 1
        assert agent.incar.release_time > 1

    model.time += 1
    model.step(tmpdir)
    model.reset_trackers()

    for agent in should_release_t2:
        assert not agent.incar.active

    # agents should not hiv convert during incar
    for agent in time_1_incar_hiv_neg:
        assert not agent.hiv.active


@pytest.mark.integration_stochastic
def test_assort_mix(params_integration, tmpdir):
    """
    Do vastly different assorting rules result in different networks
    """
    path_a = tmpdir.mkdir("a")
    path_a.mkdir("network")
    path_b = tmpdir.mkdir("b")
    path_b.mkdir("network")
    path_c = tmpdir.mkdir("c")
    path_c.mkdir("network")

    params_integration.features.assort_mix = True
    params_integration.assort_mix = ObjMap(
        {
            "same_race": {
                "attribute": "race",
                "partner_attribute": "__agent__",
                "bond_types": [],
                "agent_value": "__any__",
                "partner_values": {"__same__": 0.9, "__other__": 0.1},
            }
        }
    )

    model_a = TITAN(params_integration)
    model_a.run(path_a)

    params_integration.assort_mix = ObjMap(
        {
            "cross_race": {
                "attribute": "race",
                "partner_attribute": "__agent__",
                "bond_types": [],
                "agent_value": "__any__",
                "partner_values": {"__same__": 0.5, "__other__": 0.5},
            }
        }
    )

    model_b = TITAN(params_integration)
    model_b.run(path_b)

    params_integration.features.assort_mix = False
    model_c = TITAN(params_integration)
    model_c.run(path_c)

    model_a_rels = model_a.pop.relationships
    a_same_race = sum([1 for r in model_a_rels if r.agent1.race == r.agent2.race])
    a_diff_race = sum([1 for r in model_a_rels if r.agent1.race != r.agent2.race])

    model_b_rels = model_b.pop.relationships
    b_same_race = sum([1 for r in model_b_rels if r.agent1.race == r.agent2.race])
    b_diff_race = sum([1 for r in model_b_rels if r.agent1.race != r.agent2.race])

    model_c_rels = model_c.pop.relationships
    c_same_race = sum([1 for r in model_c_rels if r.agent1.race == r.agent2.race])
    c_diff_race = sum([1 for r in model_c_rels if r.agent1.race != r.agent2.race])

    assert a_same_race > a_diff_race
    assert a_same_race > b_same_race
    assert b_diff_race > a_diff_race
    assert math.isclose(b_same_race, c_same_race, abs_tol=50)

    # close to proportion expected
    assert math.isclose(0.9, a_same_race / (a_same_race + a_diff_race), abs_tol=0.05)
    assert math.isclose(0.5, b_same_race / (b_same_race + b_diff_race), abs_tol=0.1)
    assert math.isclose(0.5, c_same_race / (c_same_race + c_diff_race), abs_tol=0.1)


@pytest.mark.integration_stochastic
def test_treatment_cascade(params_integration, tmpdir):
    """
    Does an increase in HAART result in less HIV, AIDS and death
    """
    path_a = tmpdir.mkdir("a")
    path_a.mkdir("network")
    path_b = tmpdir.mkdir("b")
    path_b.mkdir("network")

    params_integration.features.die_and_replace = True
    params_integration.model.num_pop = 1000
    params_integration.model.time.num_steps = 20

    def run_get_death(model, path):
        deaths = set()
        hiv_start = set(a.id for a in model.pop.all_agents if a.hiv.active)
        aids_start = set(a.id for a in model.pop.all_agents if a.hiv.aids)
        hiv_only_start = hiv_start - aids_start
        unique_hiv = set()
        unique_aids = set()
        while model.time < model.params.model.time.num_steps:
            model.time += 1
            model.step(path)
            deaths |= set(a.id for a in model.deaths)
            unique_hiv |= set(
                a.id
                for a in model.pop.all_agents
                if a.hiv.active and a.hiv.time == model.time
            )
            unique_aids |= set(a.id for a in model.pop.all_agents if a.hiv.aids)
            model.reset_trackers()

        return (
            len(unique_hiv - unique_aids - hiv_only_start),
            len(unique_aids - aids_start),
            len(deaths & (unique_hiv | hiv_start)),
        )

    model_a = TITAN(params_integration)  # low haart
    new_hiv_a, new_aids_a, hiv_deaths_a = run_get_death(model_a, path_a)

    for race in params_integration.classes.races:
        for drug_type in params_integration.classes.drug_types:
            haart_params = (
                params_integration.demographics[race]
                .sex_type.MSM.drug_type[drug_type]
                .haart
            )
            haart_params.init = 0.9
            haart_params.enroll.rule.prob = 0.9
            haart_params.adherence.init = 0.9
            haart_params.adherence.prob = 0.9

    model_b = TITAN(params_integration)  # high haart
    new_hiv_b, new_aids_b, hiv_deaths_b = run_get_death(model_b, path_b)

    assert hiv_deaths_b <= hiv_deaths_a, "HIV deaths not down"
    assert new_aids_b <= new_aids_a, "new AIDS  not down"
    assert new_hiv_b < new_hiv_a, "new HIV not down"
