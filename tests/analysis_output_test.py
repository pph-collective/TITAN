import pytest

import uuid
import csv
import os
import shutil
import networkx as nx

from titan.analysis_output import *
from titan import agent
from titan.network_graph_tools import NetworkClass
from titan.params_parse import create_params


@pytest.fixture
def params(tmpdir):
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )
    return create_params(None, param_file, tmpdir)


@pytest.fixture
def stats(params):
    a = agent.Agent("MSM", 20, "BLACK", "Inj")
    a.hiv = True
    a.aids = True
    a.hiv_dx = True
    a.haart = True
    a.sne = True
    a.prep = True
    a.intervention_ever = True
    a.high_risk = True
    a.high_risk_ever = True
    a.incar = True
    a.incar_ever = True
    a.prep_reason = ["PWID", "MSMW", "HIV test"]

    p = agent.Agent("MSM", 20, "BLACK", "Inj")
    rel = agent.Relationship(a, p, 12, rel_type="sexOnly")

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
        agent_set,
        rel_list,
        agent_set,
        agent_set,
        agent_list,
        params,
    )
    return stats


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
    assert stats["ALL"]["ALL"]["numRels"] == 1
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


def test_deathReport(stats, params, tmpdir):
    run_id = uuid.uuid4()

    deathReport(run_id, 0, 1, 2, 3, stats, params, tmpdir)

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


def test_incarReport(stats, params, tmpdir):
    run_id = uuid.uuid4()

    incarReport(run_id, 0, 1, 2, 3, stats, params, tmpdir)

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


def test_newlyhighriskReport(stats, params, tmpdir):
    run_id = uuid.uuid4()

    newlyhighriskReport(run_id, 0, 1, 2, 3, stats, params, tmpdir)

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


def test_prepReport(stats, params, tmpdir):
    run_id = uuid.uuid4()

    prepReport(run_id, 0, 1, 2, 3, stats, params, tmpdir)

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


def test_basicReport(stats, params, tmpdir):
    run_id = uuid.uuid4()

    basicReport(run_id, 0, 1, 2, 3, stats, params, tmpdir)

    result_file = os.path.join(tmpdir, "basicReport_MSM_BLACK.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["rseed"] == "1"
            assert row["pseed"] == "2"
            assert row["nseed"] == "3"
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
            assert row["nseed"] == "3"
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
            assert row["nseed"] == "3"
            assert row["Total"] == "1"
            assert row["HIV"] == "1"
            assert row["PrEP"] == "1"
            assert row["Deaths"] == "1"


def test_print_components(stats, params, tmpdir):
    run_id = uuid.uuid4()

    params.model.num_pop = 1
    net = NetworkClass(params)
    components = list(net.G.subgraph(c).copy() for c in nx.connected_components(net.G))

    print_components(run_id, 0, 1, 2, 3, components, tmpdir)

    result_file = os.path.join(tmpdir, "componentReport_ALL.txt")
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["runseed"] == "1"
            assert row["compID"] == "0"
            assert row["totalN"] == "1"
