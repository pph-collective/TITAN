import pytest

from titan.network_graph_tools import *
from titan import params
from titan import agent

import os
import shutil

n_pop = 10


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
    def _make_agent(SO="MSM", age=30, race="BLACK", DU="NDU"):
        return agent.Agent(SO, age, race, DU)

    return _make_agent


def test_network_init_scale_free():
    """Test if all non-IDU,ND,NIDU agents are in the population"""
    net = NetworkClass(N=n_pop)
    assert n_pop == net.All_agentSet.num_members()

    for agent in net.All_agentSet.get_agents():
        assert agent in net.G.nodes()

    for agent in net.All_agentSet.get_agents():
        assert agent._DU in ["IDU", "NIDU", "NDU"]
        assert agent._SO in params.agentSexTypes

    assert net.get_Graph() == net.G


def test_network_init_max_k():
    """Test if all non-IDU,ND,NIDU agents are in the population"""
    net = NetworkClass(N=n_pop, network_type="max_k_comp_size")
    assert n_pop == net.All_agentSet.num_members()

    for agent in net.All_agentSet.get_agents():
        assert agent in net.G.nodes()

    for agent in net.All_agentSet.get_agents():
        assert agent._DU in ["IDU", "NIDU", "NDU"]
        assert agent._SO in params.agentSexTypes


def test_population_consistency_DU():
    """Test if Drug users add up"""
    net = NetworkClass(N=n_pop)
    check_sum_DU = (
        net.DU_IDU_agentSet.num_members()
        + net.DU_NIDU_agentSet.num_members()
        + net.DU_NDU_agentSet.num_members()
    )

    assert net.drugUse_agentSet.num_members() == check_sum_DU
    assert net.PopulationSize == check_sum_DU


def test_population_consistency_HIV():
    """Test HIV consistency"""
    net = NetworkClass(N=n_pop)
    for agent in net.All_agentSet.get_agents():
        if agent._HIV_bool:
            assert agent in net.HIV_agentSet.get_agents()

    for agent in net.HIV_agentSet.get_agents():
        assert agent._HIV_bool


def test_write_G_edgelist(setup_results_dir):
    path = "results/network/Edgelist_t0.txt"
    net = NetworkClass(N=n_pop)

    net.write_G_edgelist(path)

    count = len(open(path).readlines())

    assert count == len(net.Relationships)


def test_write_network_stats(setup_results_dir):
    path = "results/network/networkStats.txt"
    net = NetworkClass(N=n_pop)

    net.write_network_stats(path)

    asserted = False
    with open(path, "r") as f:
        for line in f:
            if "Number of nodes:" in line:
                assert int(line.split(" ")[-1]) == n_pop
                asserted = True

    # make sure we tested something was tested
    assert asserted


def test_create_graph_from_agents(make_agent):
    a = make_agent()
    b = make_agent()

    s = agent.Agent_set("test")

    s.add_agent(a)
    s.add_agent(b)

    net = NetworkClass(N=n_pop)

    assert net.G.number_of_nodes() == n_pop

    net.create_graph_from_agents(s)

    assert net.G.number_of_nodes() == n_pop + 2


def test_get_network_color():
    net = NetworkClass(N=n_pop)

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
    assert "y" in colors  # pop includes WHITE

    colors = net.get_network_color("MSW")
    assert len(colors) == n_pop
    assert "g" in colors  # pop includes WHITE
