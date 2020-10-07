import pytest

import csv
import os
import shutil
import networkx as nx
import nanoid

from titan.output import *
from titan import agent, features


@pytest.fixture
def stats(params, world_location):
    a = agent.Agent("MSM", 20, "black", "Inj", world_location)
    a.hiv = True
    a.hiv_time = 1
    a.aids = True
    a.hiv_dx = True
    a.haart.active = True
    a.syringe_services.active = True
    a.prep.active = True
    a.prep.init_class(params)
    a.prep.add_agent(a)  # add to new agents
    a.random_trial.active = True
    a.random_trial.treated = True
    a.high_risk.active = True
    a.high_risk.ever = True
    a.high_risk.add_agent(a, new_agent=True)  # add to new agents
    a.incar.active = True
    a.incar.add_agent(a)  # adds to new releases

    agent_set = agent.AgentSet("test")
    agent_set.add_agent(a)
    agent_list = [a]
    feat_list = [feature for feature in features.BaseFeature.__subclasses__()]
    stats = get_stats(agent_set, agent_set, agent_list, params, feat_list)
    return stats


@pytest.mark.unit
def test_get_stats(stats):
    assert stats["world"]["black"]["MSM"]["agents"] == 1
    assert stats["world"]["white"]["MSM"]["agents"] == 0
    assert stats["world"]["black"]["MSM"]["incar"] == 1
    assert stats["world"]["black"]["MSM"]["incar_hiv"] == 1
    assert stats["world"]["black"]["MSM"]["new_release"] == 1
    assert stats["world"]["black"]["MSM"]["new_release_hiv"] == 1
    assert stats["world"]["black"]["MSM"]["hiv_new"] == 1
    assert stats["world"]["black"]["MSM"]["hiv_new_high_risk_ever"] == 1
    assert stats["world"]["black"]["MSM"]["hiv_new_high_risk"] == 1
    assert stats["world"]["black"]["MSM"]["prep"] == 1
    assert stats["world"]["black"]["MSM"]["prep_new"] == 1
    assert stats["world"]["black"]["MSM"]["dx_new"] == 1
    assert stats["world"]["black"]["MSM"]["high_risk_new"] == 1
    assert stats["world"]["black"]["MSM"]["high_risk_new_hiv"] == 1
    assert stats["world"]["black"]["MSM"]["high_risk_new_aids"] == 1
    assert stats["world"]["black"]["MSM"]["high_risk_new_dx"] == 1
    assert stats["world"]["black"]["MSM"]["high_risk_new_haart"] == 1
    assert stats["world"]["black"]["MSM"]["hiv"] == 1
    assert stats["world"]["black"]["MSM"]["aids"] == 1
    assert stats["world"]["black"]["MSM"]["dx"] == 1
    assert stats["world"]["black"]["MSM"]["haart"] == 1
    assert stats["world"]["black"]["MSM"]["deaths"] == 1
    assert stats["world"]["black"]["MSM"]["deaths_hiv"] == 1


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
            if row["race"] == "black" and row["sex_type"] == "MSM":
                assert row["agents"] == "1"
                assert row["hiv"] == "1"
                assert row["prep"] == "1"
                assert row["deaths"] == "1"
            else:
                assert row["agents"] == "0"
                assert row["hiv"] == "0"
                assert row["prep"] == "0"
                assert row["deaths"] == "0"


@pytest.mark.unit
def test_print_components(stats, params, make_population, tmpdir):
    run_id = nanoid.generate(size=8)

    net = make_population(n=1)
    components = net.connected_components()

    print_components(run_id, 0, 1, 2, components, tmpdir, params.classes.races)

    result_file = os.path.join(tmpdir, f"{run_id}_componentReport_ALL.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["runseed"] == "1"
            assert row["compID"] == "0"
            assert row["totalN"] == "1"
