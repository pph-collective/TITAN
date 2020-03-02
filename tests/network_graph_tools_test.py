import pytest

from titan.network_graph_tools import NetworkGraphUtils
from titan.agent import Agent
from titan.population_network import PopulationClass

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
        return Agent(SO, age, race, DU)

    return _make_agent


def test_write_G_edgelist(setup_results_dir):
    path = "results/network/Edgelist_t0.txt"
    net = PopulationClass(n=n_pop, enable_nx_graph=True)
    net_util = NetworkGraphUtils(net.nx_graph)
    net_util.write_G_edgelist(path)

    count = len(open(path).readlines())

    assert count == len(net.Relationships)


def test_write_network_stats(setup_results_dir):
    path = "results/network/networkStats.txt"
    net = PopulationClass(n=n_pop, enable_nx_graph=True)

    net_util = NetworkGraphUtils(net.nx_graph)
    net_util.write_network_stats(path)

    asserted = False
    with open(path, "r") as f:
        for line in f:
            if "Number of nodes:" in line:
                assert int(line.split(" ")[-1]) == n_pop
                asserted = True

    # make sure we tested something was tested
    assert asserted


def test_get_network_color():
    net = PopulationClass(n=n_pop, enable_nx_graph=True)

    net_util = NetworkGraphUtils(net.nx_graph)
    colors = net_util.get_network_color("SO")
    assert len(colors) == n_pop
    assert "r" in colors  # pop includes MSM

    colors = net_util.get_network_color("DU")
    assert len(colors) == n_pop
    assert "g" in colors  # pop includes ND

    colors = net_util.get_network_color("Tested")
    assert len(colors) == n_pop
    assert "purple" in colors  # pop includes not HIV/tested/HAART

    colors = net_util.get_network_color("Trtmt")
    assert len(colors) == n_pop
    assert "gray" in colors  # pop includes not HIV

    colors = net_util.get_network_color("HIV")
    assert len(colors) == n_pop
    assert "g" in colors  # pop includes not HIV

    colors = net_util.get_network_color("HR")
    assert len(colors) == n_pop
    assert "g" in colors  # pop includes not HR

    colors = net_util.get_network_color("Race")
    assert len(colors) == n_pop
    assert "y" in colors  # pop includes WHITE

    colors = net_util.get_network_color("MSW")
    assert len(colors) == n_pop
    assert "g" in colors  # pop includes WHITE
