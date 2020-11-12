import pytest
import networkx as nx

import os
import subprocess
import csv
import math
from copy import deepcopy
from glob import glob

from titan.parse_params import ObjMap, create_params
from titan.model import HIVModel

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
        return HIVModel(params_integration)

    return _make_model_integration


@pytest.mark.integration_deterministic
def test_model_runs():
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )

    subprocess.check_call([f, f"-p {param_file}"])
    assert True


@pytest.mark.integration_deterministic
def test_model_reproducible(tmpdir):
    path_a = tmpdir.mkdir("result_a")
    path_b = tmpdir.mkdir("result_b")
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic_seeded.yml"
    )

    subprocess.check_call([f, f"-p {param_file}", f"-o {path_a}"])
    subprocess.check_call([f, f"-p {param_file}", f"-o {path_b}"])

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
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )

    subprocess.check_call([f, "-p", param_file, "-o", path_a, "--savepop"])

    saved_pop_path = glob(os.path.join(path_a, "pop", "*_pop.tar.gz"))[0]

    path_b = tmpdir.mkdir("b")
    subprocess.check_call(
        [f, "-p", param_file, "-o", path_b, "--poppath", saved_pop_path]
    )

    assert True


@pytest.mark.integration_deterministic
def test_model_settings_run(tmpdir):
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "integration_base.yml"
    )

    for item in os.listdir(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "settings")
    ):
        if "__" not in item and item != "base":
            path = tmpdir.mkdir(item)
            print(f"-----------Starting run for {item}-----------")
            subprocess.check_call(
                [f, f"-p {param_file}", f"-o {path}", f"-S {item}", "-e"]
            )
            assert True


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
        model_a.params.demographics.black.MSM.num_partners[bond].vars[1].value *= 10
        model_a.params.demographics.black.PWID.num_partners[bond].vars[1].value *= 10
    model_a.params.model.seed.run = model_a.run_seed
    model_a.params.model.seed.ppl = model_a.pop.pop_seed

    model_b = HIVModel(model_a.params)
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
    model_a.params.prep.target = 0.9
    model_a.params.model.seed.run = model_a.run_seed
    model_a.params.model.seed.ppl = model_a.pop.pop_seed

    model_b = HIVModel(model_a.params)
    model_b.run(path_b)
    print(
        f"model b prep world target: {model_b.pop.geography.locations['world'].params.prep.target}"
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
    assert t10_diff > t0_diff
    assert res_a["10"]["prep"] < res_b["10"]["prep"]


@pytest.mark.integration_stochastic
def test_syringe_services(make_model_integration, tmpdir):
    """
    If we use syringe services, does the incidence of hiv decrease?
    """
    model_a = make_model_integration()

    path_a = tmpdir.mkdir("a")
    path_a.mkdir("network")
    path_b = tmpdir.mkdir("b")
    path_b.mkdir("network")

    # run with default bins (0-9)
    model_a.run(path_a)

    # change the coverage upward for creating model b, use same seeds
    model_a.params.features.syringe_services = True
    model_a.params.model.seed.run = model_a.run_seed
    model_a.params.model.seed.ppl = model_a.pop.pop_seed

    model_b = HIVModel(model_a.params)

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
def test_incar(make_model_integration, tmpdir):
    model = make_model_integration()

    # turn on incar - initi is set to 0, so for these purposes, just run time
    model.params.features.incar = True

    tmpdir.mkdir("network")

    # make sure we're starting from a clean slate
    init_incar = sum([1 for agent in model.pop.all_agents if agent.incar])
    assert init_incar == 0

    model.time = 1
    model.step(tmpdir)
    model.reset_trackers()

    time_1_incars = [agent for agent in model.pop.all_agents if agent.incar]
    time_1_incar_hiv_neg = [agent for agent in time_1_incars if not agent.hiv]
    should_release_t2 = [agent for agent in time_1_incars if agent.incar_time == 1]

    assert len(time_1_incars) > 0
    assert len(time_1_incar_hiv_neg) > 0
    assert len(should_release_t2) > 0

    for agent in time_1_incars:
        assert agent.incar_time > 0

    model.time += 1
    model.step(tmpdir)

    for agent in should_release_t2:
        assert agent in model.new_incar_release
        assert not agent.incar
        assert agent.incar_time == 0

    model.reset_trackers()

    # agents should not hiv convert during incar
    for agent in time_1_incar_hiv_neg:
        assert not agent.hiv


@pytest.mark.integration_deterministic
def test_incar(make_model_integration, tmpdir):
    model = make_model_integration()

    # turn on partner tracing - just run time affects
    model.params.features.partner_tracing = True

    tmpdir.mkdir("network")
