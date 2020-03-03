import pytest

from titan.network import *
from titan.parse_params import create_params
from titan import agent

import os
import shutil

n_pop = 100


@pytest.fixture
def setup_results_dir():
    outfile_dir = os.path.join(os.getcwd(), "results", "network")
    if os.path.isdir(outfile_dir):
        shutil.rmtree(outfile_dir)
    os.makedirs(outfile_dir)
    yield outfile_dir
    shutil.rmtree(outfile_dir)


@pytest.fixture
def make_agent():
    def _make_agent(SO="MSM", age=30, race="BLACK", DU="None"):
        return agent.Agent(SO, age, race, DU)

    return _make_agent


@pytest.fixture
def params(tmpdir):
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )
    return create_params(None, param_file, tmpdir)


def test_network_init_scale_free(params):
    """Test if all Inj,NonInj,None drug use agents are in the population"""
    net = Network(params)
    assert n_pop == net.all_agents.num_members()

    for agent in net.all_agents:
        assert agent in net.G.nodes()

    for agent in net.all_agents:
        assert agent.drug_use in ["Inj", "NonInj", "None"]
        assert agent.so in params.classes.sex_types

    assert net.G == net.G


def test_network_init_max_k(params):
    """Test if all Inj,NonInj,None drug use agents are in the population"""
    params.model.network.type = "max_k_comp_size"
    net = Network(params)
    assert n_pop == net.all_agents.num_members()

    for agent in net.all_agents:
        assert agent in net.G.nodes()

    for agent in net.all_agents:
        assert agent.drug_use in ["Inj", "NonInj", "None"]
        assert agent.so in params.classes.sex_types


def test_population_consistency_DU(params):
    """Test if Drug users add up"""
    net = Network(params)
    check_sum_DU = (
        net.drug_use_inj_agents.num_members()
        + net.drug_use_noninj_agents.num_members()
        + net.drug_use_none_agents.num_members()
    )

    assert net.drug_use_agents.num_members() == check_sum_DU
    assert params.model.num_pop == check_sum_DU


def test_population_consistency_HIV(params):
    """Test HIV consistency"""
    net = Network(params)
    for agent in net.all_agents:
        if agent.hiv:
            assert agent in net.hiv_agents

    for agent in net.hiv_agents:
        assert agent.hiv


def test_write_graph_edgelist(setup_results_dir, params):
    path = "results/network/Edgelist_t0.txt"
    net = Network(params)

    net.write_graph_edgelist(path)

    count = len(open(path).readlines())

    assert count == len(net.relationships)


def test_write_network_stats(setup_results_dir, params):
    path = "results/network/networkStats.txt"
    net = Network(params)

    net.write_network_stats(path)

    asserted = False
    with open(path, "r") as f:
        for line in f:
            if "Number of nodes:" in line:
                assert int(line.split(" ")[-1]) == n_pop
                asserted = True

    # make sure we tested something was tested
    assert asserted


def test_create_graph_from_agents(make_agent, params):
    a = make_agent()
    b = make_agent()

    s = agent.AgentSet("test")

    s.add_agent(a)
    s.add_agent(b)

    net = Network(params)

    assert net.G.number_of_nodes() == n_pop

    net.create_graph_from_agents(s)

    assert net.G.number_of_nodes() == n_pop + 2


def test_get_network_color(params):
    net = Network(params)

    colors = net.get_network_color("SO")
    assert len(colors) == n_pop
    assert "r" in colors  # pop includes MSM

    colors = net.get_network_color("DU")
    assert len(colors) == n_pop
    assert "g" in colors  # pop includes ND

    colors = net.get_network_color("Tested")
    assert len(colors) == n_pop
    assert "purple" in colors  # pop includes not HIV/tested/HAART

    colors = net.get_network_color("Trtmt")
    assert len(colors) == n_pop
    assert "gray" in colors  # pop includes not HIV

    colors = net.get_network_color("HIV")
    assert len(colors) == n_pop
    assert "g" in colors  # pop includes not HIV

    colors = net.get_network_color("HR")
    assert len(colors) == n_pop
    assert "g" in colors  # pop includes not HR

    colors = net.get_network_color("Race")
    assert len(colors) == n_pop
    assert "g" in colors  # pop includes BLACK

    colors = net.get_network_color("MSW")
    assert len(colors) == n_pop
    assert "y" in colors  # pop includes BLACK
