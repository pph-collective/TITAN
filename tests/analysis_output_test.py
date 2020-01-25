import pytest

import uuid
import csv
import os
import shutil
import networkx as nx

from titan.analysis_output import *
from titan import agent
from titan.network_graph_tools import NetworkClass


@pytest.fixture
def setup_results_dir():
    outfile_dir = os.path.join(os.getcwd(), "results")
    if os.path.isdir(outfile_dir):
        shutil.rmtree(outfile_dir)
    os.mkdir(outfile_dir)
    yield outfile_dir
    shutil.rmtree(outfile_dir)


@pytest.fixture
def stats():
    a = agent.Agent("MSM", 20, "BLACK", "IDU")
    a._HIV_bool = True
    a._AIDS_bool = True
    a._diagnosed = True
    a._HAART_bool = True
    a._SNE_bool = True
    a._PrEP_bool = True
    a._treatment_bool = True
    a._highrisk_bool = True
    a._everhighrisk_bool = True
    a._incar_bool = True
    a._ever_incar_bool = True
    a._PrEP_reason = ["IDU", "MSMW", "HIV test"]

    p = agent.Agent("MSM", 20, "BLACK", "IDU")
    rel = agent.Relationship(a, p, 12, rel_type="sexOnly")

    agent_set = agent.Agent_set("test")
    DU_set = agent.Agent_set("DU", agent_set)
    IDU_set = agent.Agent_set("IDU", DU_set)
    IDU_set.add_agent(a)
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
    assert stats["BLACK"]["MSM"]["newHR_diagnosed"] == 1
    assert stats["BLACK"]["MSM"]["newHR_ART"] == 1
    assert stats["BLACK"]["MSM"]["numHIV"] == 1
    assert stats["BLACK"]["MSM"]["numAIDS"] == 1
    assert stats["BLACK"]["MSM"]["numTested"] == 1
    assert stats["BLACK"]["MSM"]["numART"] == 1
    assert stats["BLACK"]["IDU"]["numAgents"] == 1
    assert stats["BLACK"]["IDU"]["numHIV"] == 1
    assert stats["BLACK"]["IDU"]["numAIDS"] == 1
    assert stats["BLACK"]["IDU"]["numTested"] == 1
    assert stats["BLACK"]["IDU"]["numART"] == 1
    assert stats["BLACK"]["MSM"]["deaths"] == 1
    assert stats["BLACK"]["MSM"]["deaths_HIV"] == 1
    assert stats["ALL"]["ALL"]["numRels"] == 1
    assert stats["ALL"]["ALL"]["numAgents"] == 1
    assert stats["BLACK"]["ALL"]["numAgents"] == 1
    assert stats["ALL"]["MSM"]["numAgents"] == 1
    assert stats["ALL"]["HM"]["numAgents"] == 0
    assert stats["ALL"]["IDU"]["numAgents"] == 1
    assert stats["BLACK"]["IDU"]["numAgents"] == 1
    assert stats["BLACK"]["IDU"]["numART"] == 1
    assert stats["ALL"]["ALL"]["numART"] == 1
    assert stats["ALL"]["IDU"]["numART"] == 1
    assert stats["BLACK"]["ALL"]["numART"] == 1


def test_deathReport(stats, setup_results_dir):
    run_id = uuid.uuid4()

    deathReport(run_id, 0, 1, 2, 3, stats)

    result_file = "results/DeathReport.txt"
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


def test_incarReport(stats, setup_results_dir):
    run_id = uuid.uuid4()

    incarReport(run_id, 0, 1, 2, 3, stats)

    result_file = "results/IncarReport.txt"
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


def test_newlyhighriskReport(stats, setup_results_dir):
    run_id = uuid.uuid4()

    newlyhighriskReport(run_id, 0, 1, 2, 3, stats)

    result_file = "results/newlyHR_Report.txt"
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


def test_prepReport(stats, setup_results_dir):
    run_id = uuid.uuid4()

    prepReport(run_id, 0, 1, 2, 3, stats)

    result_file = "results/PrEPReport.txt"
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["seed"] == "1"
            assert row["NewEnroll"] == "1"
            assert row["IDUpartner"] == "1"
            assert row["TestedPartner"] == "1"
            assert row["MSMWpartner"] == "1"


def test_basicReport(stats, setup_results_dir):
    run_id = uuid.uuid4()

    basicReport(run_id, 0, 1, 2, 3, stats)

    result_file = "results/basicReport_MSM_BLACK.txt"
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

    result_file = "results/basicReport_HM_WHITE.txt"
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

    result_file = "results/basicReport_ALL_ALL.txt"
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


def test_print_components(stats, setup_results_dir):
    run_id = uuid.uuid4()

    net = NetworkClass(N=1)
    G = net.get_Graph()
    components = list(G.subgraph(c).copy() for c in nx.connected_components(G))

    print_components(run_id, 0, 1, 2, 3, components)

    result_file = "results/componentReport_ALL.txt"
    assert os.path.isfile(result_file)
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            assert row["t"] == "0"
            assert row["run_id"] == str(run_id)
            assert row["runseed"] == "1"
            assert row["compID"] == "0"
            assert row["totalN"] == "1"
