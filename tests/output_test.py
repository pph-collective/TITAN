import pytest

import csv
import os
import shutil
import networkx as nx
import nanoid

from titan.output import *
from titan import agent


@pytest.fixture
def stats(params):
    a = agent.Agent("MSM", 20, "BLACK", "Inj")
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
    a.prep_reason = ["PWID", "MSMW", "HIV test"]

    p = agent.Agent("MSM", 20, "BLACK", "Inj")
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
    assert stats["WHITE"]["ALL"]["numAgents"] == 0
    assert stats["BLACK"]["MSM"]["numAgents"] == 1
    assert stats["BLACK"]["ALL"]["numAgents"] == 1
    assert stats["BLACK"]["MSM"]["incar"] == 1
    assert stats["BLACK"]["MSM"]["incarHIV"] == 1
    assert stats["BLACK"]["MSM"]["newRelease"] == 1
    assert stats["BLACK"]["MSM"]["newReleaseHIV"] == 1
    assert stats["BLACK"]["MSM"]["inf_newInf"] == 1
    assert stats["BLACK"]["MSM"]["inf_HRever"] == 1
    assert stats["BLACK"]["MSM"]["inf_HR6m"] == 1
    assert stats["BLACK"]["MSM"]["numPrEP"] == 1
    assert stats["BLACK"]["MSM"]["iduPartPrep"] == 1
    assert stats["BLACK"]["MSM"]["msmwPartPrep"] == 1
    assert stats["BLACK"]["MSM"]["testedPartPrep"] == 1
    assert stats["BLACK"]["MSM"]["newNumPrEP"] == 1
    assert stats["BLACK"]["MSM"]["newlyTested"] == 1
    assert stats["BLACK"]["MSM"]["newHR"] == 1
    assert stats["ALL"]["MSM"]["newHR"] == 1
    assert stats["BLACK"]["MSM"]["newHR_HIV"] == 1
    assert stats["BLACK"]["MSM"]["newHR_AIDS"] == 1
    assert stats["BLACK"]["MSM"]["newHR_tested"] == 1
    assert stats["BLACK"]["MSM"]["newHR_ART"] == 1
    assert stats["BLACK"]["MSM"]["numHIV"] == 1
    assert stats["BLACK"]["MSM"]["numAIDS"] == 1
    assert stats["BLACK"]["MSM"]["numTested"] == 1
    assert stats["BLACK"]["MSM"]["numART"] == 1
    assert stats["BLACK"]["PWID"]["numAgents"] == 1
    assert stats["BLACK"]["PWID"]["numHIV"] == 1
    assert stats["BLACK"]["PWID"]["numAIDS"] == 1
    assert stats["BLACK"]["PWID"]["numTested"] == 1
    assert stats["BLACK"]["PWID"]["numART"] == 1
    assert stats["BLACK"]["MSM"]["deaths"] == 1
    assert stats["BLACK"]["MSM"]["deaths_HIV"] == 1
    assert stats["ALL"]["ALL"]["numAgents"] == 1
    assert stats["BLACK"]["ALL"]["numAgents"] == 1
    assert stats["ALL"]["MSM"]["numAgents"] == 1
    assert stats["ALL"]["HM"]["numAgents"] == 0
    assert stats["ALL"]["PWID"]["numAgents"] == 1
    assert stats["BLACK"]["PWID"]["numAgents"] == 1
    assert stats["BLACK"]["PWID"]["numART"] == 1
    assert stats["ALL"]["ALL"]["numART"] == 1
    assert stats["ALL"]["PWID"]["numART"] == 1
    assert stats["BLACK"]["ALL"]["numART"] == 1


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
            assert row["seed"] == "1"
            assert row["tot_MSM"] == "1"
            assert row["tot_HM"] == "0"
            assert row["HIV_MSM"] == "1"
            assert row["HIV_HM"] == "0"


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
            assert row["seed"] == "1"
            assert row["ALL_MSM_tot"] == "1"
            assert row["BLACK_MSM_HIV"] == "1"
            assert row["ALL_HM_rlsd"] == "0"
            assert row["WHITE_MSM_rlsdHIV"] == "0"


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
            assert row["seed"] == "1"
            assert row["newHR_HM"] == "0"
            assert row["newHR_MSM"] == "1"
            assert row["newHR_Tested_MSM"] == "1"
            assert row["newHR_ART_MSM"] == "1"


@pytest.mark.unit
def test_prepReport(stats, params, tmpdir):
    run_id = nanoid.generate(size=8)

    prepReport(run_id, 0, 1, 2, stats, params, tmpdir)

    result_file = os.path.join(tmpdir, "PrEPReport.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["seed"] == "1"
            assert row["NewEnroll"] == "1"
            assert row["PWIDpartner"] == "1"
            assert row["TestedPartner"] == "1"
            assert row["MSMWpartner"] == "1"


@pytest.mark.unit
def test_basicReport(stats, params, tmpdir):
    run_id = nanoid.generate(size=8)

    basicReport(run_id, 0, 1, 2, stats, params, tmpdir)

    result_file = os.path.join(tmpdir, "basicReport_MSM_BLACK.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["rseed"] == "1"
            assert row["pseed"] == "2"
            assert row["Total"] == "1"
            assert row["HIV"] == "1"
            assert row["PrEP"] == "1"
            assert row["Deaths"] == "1"

    result_file = os.path.join(tmpdir, "basicReport_HM_WHITE.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["rseed"] == "1"
            assert row["pseed"] == "2"
            assert row["Total"] == "0"
            assert row["HIV"] == "0"
            assert row["PrEP"] == "0"
            assert row["Deaths"] == "0"

    result_file = os.path.join(tmpdir, "basicReport_ALL_ALL.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["rseed"] == "1"
            assert row["pseed"] == "2"
            assert row["Total"] == "1"
            assert row["HIV"] == "1"
            assert row["PrEP"] == "1"
            assert row["Deaths"] == "1"


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
