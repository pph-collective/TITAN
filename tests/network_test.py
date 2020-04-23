import pytest

from titan.network import *
from titan.parse_params import create_params
from titan import agent
from titan.population import Population

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


def test_write_graph_edgelist(setup_results_dir, params):
    path = "results/network/Edgelist_t0.txt"
    net = Population(params)
    net_util = NetworkGraphUtils(net.graph)

    net_util.write_graph_edgelist(path)

    count = len(open(path).readlines())

    assert count == len(net.relationships)


def test_write_network_stats(setup_results_dir, params):
    path = "results/network/networkStats.txt"
    net = Population(params)
    net_util = NetworkGraphUtils(net.graph)

    net_util.write_network_stats(path)

    asserted = False
    with open(path, "r") as f:
        for line in f:
            if "Number of nodes:" in line:
                assert int(line.split(" ")[-1]) == n_pop
                asserted = True

    # make sure we tested something was tested
    assert asserted


def test_get_network_color(params):
    net = Population(params)
    net_util = NetworkGraphUtils(net.graph)

    colors = net_util.get_network_color("so")
    assert len(colors) == n_pop
    assert "b" in colors  # pop includes MSM

    colors = net_util.get_network_color("drug_use")
    assert len(colors) == n_pop
    assert "b" in colors  # pop includes ND

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

    colors = net_util.get_network_color("race")
    assert len(colors) == n_pop
    assert "b" in colors  # pop includes BLACK
