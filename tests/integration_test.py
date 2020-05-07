import pytest
import networkx as nx

import os
import subprocess
import csv
import math
from copy import deepcopy

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

    result_file_a = os.path.join(path_a, "basicReport_MSM_BLACK.txt")
    result_file_b = os.path.join(path_b, "basicReport_MSM_BLACK.txt")
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
        assert res_a[i]["Total"] == res_b[i]["Total"]
        assert res_a[i]["HIV"] == res_b[i]["HIV"]
        assert res_a[i]["PrEP"] == res_b[i]["PrEP"]
        assert res_a[i]["Deaths"] == res_b[i]["Deaths"]
        assert res_a[i]["HIV"] == res_b[i]["HIV"]


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
            subprocess.check_call([f, f"-p {param_file}", f"-o {path}", f"-S {item}"])
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

    # change the bins upward for creating model b
    bins = {
        25: {"prob": 0.3},
        30: {"prob": 0.4},
        35: {"prob": 0.3},
    }
    model_a.params.model.population.num_partners.bins = ObjMap(bins)
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
    assert (g_a_0.number_of_edges() * 2) < g_b_0.number_of_edges()
    assert (g_a_10.number_of_edges() * 2) < g_b_10.number_of_edges()


@pytest.mark.integration_stochastic
def test_prep_coverage(make_model_integration, tmpdir):
    """
    If we increase the target of prep coverage, does the incidence of hiv decrease?
    """
    model_a = make_model_integration()

    path_a = tmpdir.mkdir("a")
    path_a.mkdir("network")
    path_b = tmpdir.mkdir("b")
    path_b.mkdir("network")

    # run with default bins (0-9)
    model_a.run(path_a)

    # change the coverage upward for creating model b, use same seeds
    model_a.params.prep.target = 0.9
    model_a.params.model.seed.run = model_a.run_seed
    model_a.params.model.seed.ppl = model_a.pop.pop_seed

    model_b = HIVModel(model_a.params)
    model_b.run(path_b)

    result_file_a = os.path.join(path_a, "basicReport_ALL_ALL.txt")
    result_file_b = os.path.join(path_b, "basicReport_ALL_ALL.txt")
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

    # at start, should be similar
    assert res_a[0]["t"] == res_b[0]["t"]
    t0_hiv_a = float(res_a[0]["HIV"])
    t0_hiv_b = float(res_b[0]["HIV"])
    t0_diff = t0_hiv_a - t0_hiv_b
    assert math.isclose(t0_hiv_a, t0_hiv_b, abs_tol=15)  # within 15 agents
    assert int(res_a[0]["PrEP"]) < int(res_b[0]["PrEP"])

    # at end, should be different
    assert res_a[9]["t"] == res_b[9]["t"]
    t9_hiv_a = float(res_a[9]["HIV"])
    t9_hiv_b = float(res_b[9]["HIV"])
    t9_diff = t9_hiv_a - t9_hiv_b  # a should be higher
    assert t9_diff > t0_diff
    assert int(res_a[9]["PrEP"]) < int(res_b[9]["PrEP"])


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

    result_file_a = os.path.join(path_a, "basicReport_PWID_ALL.txt")
    result_file_b = os.path.join(path_b, "basicReport_PWID_ALL.txt")
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

    # at start, should be similar
    assert res_a[0]["t"] == res_b[0]["t"]
    t0_hiv_a = float(res_a[0]["HIV"])
    t0_hiv_b = float(res_b[0]["HIV"])
    t0_diff = t0_hiv_a - t0_hiv_b
    assert math.isclose(t0_hiv_a, t0_hiv_b, abs_tol=15)  # within 15 agents

    # at end, should be different
    assert res_a[9]["t"] == res_b[9]["t"]
    t9_hiv_a = float(res_a[9]["HIV"])
    t9_hiv_b = float(res_b[9]["HIV"])
    t9_diff = t9_hiv_a - t9_hiv_b  # a should be higher
    assert t9_diff > t0_diff


@pytest.mark.integration_deterministic
def test_static_network(make_model_integration, tmpdir):
    model = make_model_integration()
    orig_rel_ids = [rel.id for rel in model.pop.relationships]

    # turn on static network - only affects run time, so fine have false during init
    model.params.features.static_network = True

    tmpdir.mkdir("network")

    for t in range(1, 10):
        model.time = t
        model.step(tmpdir, burn=False)
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
    model.step(tmpdir, burn=False)
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
    model.step(tmpdir, burn=False)

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
