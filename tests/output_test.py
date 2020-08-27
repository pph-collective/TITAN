import pytest

import csv
import os
import shutil
import networkx as nx
import nanoid

from titan.output import *
from titan import agent


@pytest.fixture
def stats(params, world_location):
    a = agent.Agent("MSM", 20, "black", "Inj", world_location)
    a.hiv = True
    a.aids = True
    a.hiv_dx = True
    a.haart = True
    a.ssp = True
    a.prep = True
    a.intervention_ever = True
    a.high_risk = True
    a.high_risk_ever = True
    a.incar = True

    p = agent.Agent("MSM", 20, "black", "Inj", world_location)
    p.partners["Sex"] = set()
    a.partners["Sex"] = set()
    rel = agent.Relationship(a, p, 12, bond_type="Sex")

    agent_set = agent.AgentSet("test")
    DU_set = agent.AgentSet("DU", agent_set)
    PWID_set = agent.AgentSet("Inj", DU_set)
    PWID_set.add_agent(a)
    agent_list = [a]
    rel_list = [rel]
    stats = get_stats(
        agent_set,
        agent_set,
        agent_set,
        agent_set,
        agent_set,
        agent_set,
        agent_list,
        params,
    )
    return stats


@pytest.mark.unit
def test_get_stats(stats):
    assert stats["world"]["black"]["MSM"]["numAgents"] == 1
    assert stats["world"]["white"]["MSM"]["numAgents"] == 0
    assert stats["world"]["black"]["MSM"]["incar"] == 1
    assert stats["world"]["black"]["MSM"]["incarHIV"] == 1
    assert stats["world"]["black"]["MSM"]["newRelease"] == 1
    assert stats["world"]["black"]["MSM"]["newReleaseHIV"] == 1
    assert stats["world"]["black"]["MSM"]["inf_newInf"] == 1
    assert stats["world"]["black"]["MSM"]["inf_HRever"] == 1
    assert stats["world"]["black"]["MSM"]["inf_HR6m"] == 1
    assert stats["world"]["black"]["MSM"]["numPrEP"] == 1
    assert stats["world"]["black"]["MSM"]["newNumPrEP"] == 1
    assert stats["world"]["black"]["MSM"]["newlyDiagnosed"] == 1
    assert stats["world"]["black"]["MSM"]["newHR"] == 1
    assert stats["world"]["black"]["MSM"]["newHR_HIV"] == 1
    assert stats["world"]["black"]["MSM"]["newHR_AIDS"] == 1
    assert stats["world"]["black"]["MSM"]["newHR_dx"] == 1
    assert stats["world"]["black"]["MSM"]["newHR_ART"] == 1
    assert stats["world"]["black"]["MSM"]["numHIV"] == 1
    assert stats["world"]["black"]["MSM"]["numAIDS"] == 1
    assert stats["world"]["black"]["MSM"]["numDiagnosed"] == 1
    assert stats["world"]["black"]["MSM"]["numART"] == 1
    assert stats["world"]["black"]["MSM"]["deaths"] == 1
    assert stats["world"]["black"]["MSM"]["deaths_HIV"] == 1


@pytest.mark.unit
def test_deathReport(stats, params, tmpdir):
    run_id = nanoid.generate(size=8)

    deathReport(run_id, 0, 1, 2, stats, params, tmpdir)

    result_file = os.path.join(tmpdir, "DeathReport.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["rseed"] == "1"
            if row["race"] == "black" and row["sex_type"] == "MSM":
                assert row["tot"] == "1"
                assert row["HIV"] == "1"
            else:
                assert row["tot"] == "0"
                assert row["HIV"] == "0"


@pytest.mark.unit
def test_incarReport(stats, params, tmpdir):
    run_id = nanoid.generate(size=8)

    incarReport(run_id, 0, 1, 2, stats, params, tmpdir)

    result_file = os.path.join(tmpdir, "IncarReport.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["rseed"] == "1"
            if row["race"] == "black" and row["sex_type"] == "MSM":
                assert row["tot"] == "1"
                assert row["HIV"] == "1"
                assert row["rlsd"] == "1"
                assert row["rlsdHIV"] == "1"
            else:
                assert row["tot"] == "0"
                assert row["HIV"] == "0"
                assert row["rlsd"] == "0"
                assert row["rlsdHIV"] == "0"


@pytest.mark.unit
def test_newlyhighriskReport(stats, params, tmpdir):
    run_id = nanoid.generate(size=8)

    newlyhighriskReport(run_id, 0, 1, 2, stats, params, tmpdir)

    result_file = os.path.join(tmpdir, "newlyHR_Report.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["rseed"] == "1"
            if row["race"] == "black" and row["sex_type"] == "MSM":
                assert row["newHR"] == "1"
                assert row["newHR_HIV"] == "1"
                assert row["newHR_AIDS"] == "1"
                assert row["newHR_Diagnosed"] == "1"
                assert row["newHR_ART"] == "1"
            else:
                assert row["newHR"] == "0"
                assert row["newHR_HIV"] == "0"
                assert row["newHR_AIDS"] == "0"
                assert row["newHR_Diagnosed"] == "0"
                assert row["newHR_ART"] == "0"


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
                assert row["Total"] == "1"
                assert row["HIV"] == "1"
                assert row["PrEP"] == "1"
                assert row["Deaths"] == "1"
            else:
                assert row["Total"] == "0"
                assert row["HIV"] == "0"
                assert row["PrEP"] == "0"
                assert row["Deaths"] == "0"


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
