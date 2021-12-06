import pytest

import csv
import os
import shutil
import networkx as nx
import nanoid

from titan.output import *
from titan import agent, features, exposures


@pytest.fixture
def stats(params, world_location):
    cur_time = 3
    a = agent.Agent("MSM", 20, "black", "Inj", world_location)
    a.hiv.active = True
    a.hiv.time = cur_time
    a.hiv.aids = True
    a.hiv.dx = True
    a.hiv.dx_time = cur_time
    a.haart.active = True
    a.syringe_services.active = True
    a.prep.active = True
    a.prep.time = cur_time
    a.random_trial.active = True
    a.random_trial.treated = True
    a.high_risk.active = True
    a.high_risk.ever = True
    a.high_risk.time = cur_time
    a.incar.active = True
    a.incar.release_time = cur_time
    a.component = 0

    agent_set = agent.AgentSet("test")
    agent_set.add_agent(a)
    agent_exit = {"death": [a], "age_out": [a]}
    feat_list = [feature for feature in features.BaseFeature.__subclasses__()]
    expose_list = [exposure for exposure in exposures.BaseExposure.__subclasses__()]
    params.classes.components = ["-1", "0"]
    stats = get_stats(agent_set, agent_exit, params, feat_list, expose_list, cur_time)
    return stats


@pytest.mark.unit
def test_get_stats(stats):
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["agents"] == 1
    assert stats["world"]["white"]["MSM"]["0"]["Inj"]["agents"] == 0
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["incar"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["incar_hiv"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["new_release"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["new_release_hiv"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["hiv_new"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["hiv_new_high_risk_ever"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["hiv_new_high_risk"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["prep"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["prep_new"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["hiv_dx_new"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["high_risk_new"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["high_risk_new_hiv"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["high_risk_new_aids"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["high_risk_new_dx"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["high_risk_new_haart"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["hiv"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["hiv_aids"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["hiv_dx"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["haart"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["death"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["age_out"] == 1
    assert stats["world"]["black"]["MSM"]["0"]["Inj"]["deaths_hiv"] == 1


@pytest.mark.unit
def test_basicReport(stats, params, tmpdir):
    run_id = nanoid.generate(size=8)

    basicReport(run_id, 0, 1, 2, stats, params, tmpdir)

    result_file = os.path.join(tmpdir, "basicReport.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["rseed"] == "1"
            if (
                row["race"] == "black"
                and row["sex_type"] == "MSM"
                and row["component"] == "0"
                and row["drug_type"] == "Inj"
            ):
                assert row["agents"] == "1"
                assert row["hiv"] == "1"
                assert row["prep"] == "1"
                assert row["death"] == "1"
            else:
                assert row["agents"] == "0"
                assert row["hiv"] == "0"
                assert row["prep"] == "0"
                assert row["deaths"] == "0"


@pytest.mark.unit
def test_print_components(stats, params, make_population, tmpdir):
    run_id = nanoid.generate(size=8)

    pop = make_population(n=1)
    components = pop.connected_components()

    print_components(run_id, 0, 1, 2, components, tmpdir)

    result_file = os.path.join(tmpdir, f"{run_id}_componentReport_ALL.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["runseed"] == "1"
            assert row["component"] == "0"


@pytest.mark.unit
def test_write_graph_edgelist(setup_results_dir, make_population):
    id = "test"
    t = 0
    path = "results/network"
    n_pop = 100
    pop = make_population(n=n_pop)

    write_graph_edgelist(pop.graph, path, id, t)

    file_path = os.path.join(path, f"{id}_Edgelist_t{t}.txt")
    count = len(open(file_path).readlines())

    assert count == len(pop.relationships)


@pytest.mark.unit
def test_write_network_stats(setup_results_dir, make_population):
    id = "test"
    t = 0
    path = "results/network"
    n_pop = 100
    pop = make_population(n=n_pop)

    write_network_stats(pop.graph, path, id, t)

    file_path = os.path.join(path, f"{id}_NetworkStats_t{t}.txt")
    asserted = False
    with open(file_path, "r") as f:
        for line in f:
            if "Number of nodes:" in line:
                assert int(line.split(" ")[-1]) == n_pop
                asserted = True

    # make sure we tested something was tested
    assert asserted
