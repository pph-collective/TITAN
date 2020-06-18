import pytest

from titan.network import *
from titan.parse_params import create_params
from titan import agent
from titan.population import Population

import os
import shutil

n_pop = 100


@pytest.mark.unit
def test_write_graph_edgelist(setup_results_dir, params):
    id = "test"
    t = 0
    path = "results/network"
    net = Population(params)
    net_util = NetworkGraphUtils(net.graph)

    net_util.write_graph_edgelist(path, id, t)

    file_path = os.path.join(path, f"{id}_Edgelist_t{t}.txt")
    count = len(open(file_path).readlines())

    assert count == len(net.relationships)


@pytest.mark.unit
def test_write_network_stats(setup_results_dir, params):
    id = "test"
    t = 0
    path = "results/network"
    net = Population(params)
    net_util = NetworkGraphUtils(net.graph)

    net_util.write_network_stats(path, id, t)

    file_path = os.path.join(path, f"{id}_NetworkStats_t{t}.txt")
    asserted = False
    with open(file_path, "r") as f:
        for line in f:
            if "Number of nodes:" in line:
                assert int(line.split(" ")[-1]) == n_pop
                asserted = True

    # make sure we tested something was tested
    assert asserted


@pytest.mark.unit
def test_get_network_color(params):
    net = Population(params)
    net_util = NetworkGraphUtils(net.graph)

    colors = net_util.get_network_color("sex_type")
    assert len(colors) == n_pop
    assert "b" in colors  # pop includes MSM

    colors = net_util.get_network_color("drug_type")
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
    assert "b" in colors  # pop includes black
