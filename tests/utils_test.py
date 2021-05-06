import pytest

import os
import sys
import random
import subprocess

import numpy as np
import oyaml as yaml

import titan.utils as utils
from titan.parse_params import ObjMap

from conftest import FakeRandom


@pytest.mark.unit
def test_check_rand_int():
    assert utils.get_check_rand_int(1) == 1
    assert utils.get_check_rand_int(0) > 0

    with pytest.raises(ValueError):
        utils.get_check_rand_int(3.3)


@pytest.mark.unit
def test_safe_divide():
    assert utils.safe_divide(1, 0) == 0.0
    assert utils.safe_divide(1, 2) == 0.5


@pytest.mark.unit
def test_safe_random_choice():
    rand_gen = random.Random(123)
    assert utils.safe_random_choice([], rand_gen) is None

    assert utils.safe_random_choice({1, 2, 3}, rand_gen) in (1, 2, 3)
    assert utils.safe_random_choice([1, 2, 3], rand_gen) == 1

    assert utils.safe_random_choice([1, 2, 3], rand_gen, [1000, 1, 1]) == 1


@pytest.mark.unit
def test_safe_shuffle():
    rand_gen = random.Random(123)
    assert utils.safe_shuffle([], rand_gen) == []

    assert utils.safe_shuffle({1, 2, 3}, rand_gen) != [1, 2, 3]
    assert utils.safe_shuffle([1, 2, 3], rand_gen) != [1, 2, 3]


@pytest.mark.unit
def test_safe_dist():
    rand_gen = np.random.RandomState(123)

    dist_info = ObjMap(
        {"dist_type": "poisson", "vars": {1: {"value": 20, "value_type": "int"}}}
    )

    # numpy dist
    assert utils.safe_dist(dist_info, rand_gen) > 0

    # custom dist
    dist_info["dist_type"] = "set_value"
    assert utils.safe_dist(dist_info, rand_gen) == 20

    # not a real distribution
    dist_info["dist_type"] = "fake_dist"
    with pytest.raises(AttributeError) as excinfo:
        utils.safe_dist(dist_info, rand_gen)
    assert "Distribution type fake_dist not found!" in str(excinfo)

    # pert
    low = 2
    peak = 5
    high = 100
    temp = 4
    pert_info = ObjMap(
        {
            "dist_type": "pert",
            "vars": {
                1: {"value": low, "value_type": "int"},
                2: {"value": peak, "value_type": "int"},
                3: {"value": high, "value_type": "int"},
                4: {"value": temp, "value_type": "int"},
            },
        }
    )

    val = utils.safe_dist(pert_info, rand_gen)
    assert val > low
    assert val < high

    # weibull_modified (2-param weibull)
    rand_gen = FakeRandom(0.5)
    shape = 0.5
    scale = 11
    weibull_info = ObjMap(
        {
            "dist_type": "weibull_modified",
            "vars": {
                1: {"value": shape, "value_type": "float"},
                2: {"value": scale, "value_type": "float"},
            },
        }
    )
    assert utils.safe_dist(weibull_info, rand_gen) == 5.284983153100216


@pytest.mark.unit
def test_get_param_from_path(params):
    assert params.classes.sex_types.HM.cis_trans == "cis"

    param_path_pipe = "classes|sex_types|HM|cis_trans"
    path_params_pipe, last_item_pipe = utils.get_param_from_path(
        params, param_path_pipe, "|"
    )
    assert last_item_pipe == "cis_trans"
    assert "cis_trans" in path_params_pipe
    assert path_params_pipe["cis_trans"] == "cis"

    param_path_hash = "classes#sex_types#HM#cis_trans"
    path_params_hash, last_item_hash = utils.get_param_from_path(
        params, param_path_hash, "#"
    )
    assert last_item_hash == "cis_trans"
    assert "cis_trans" in path_params_hash
    assert path_params_hash["cis_trans"] == "cis"

    assert params.partnership.sex.haart_scaling.HM.non_adherent == 1.0
    param_path_haart = "partnership|sex|haart_scaling|HM|non_adherent"
    path_params_haart, last_item_haart = utils.get_param_from_path(
        params, param_path_haart, "|"
    )
    assert last_item_haart == "non_adherent"
    assert "non_adherent" in path_params_haart
    assert path_params_haart["non_adherent"] == 1.0


@pytest.mark.unit
def test_scale_param(params):
    assert params.classes.locations.world.ppl == 1.0
    param_path_pipe = "classes|locations|world|ppl"

    utils.scale_param(params, param_path_pipe, 2)

    assert params.classes.locations.world.ppl == 2.0

    param_path_hash = "classes#locations#world#ppl"

    utils.scale_param(params, param_path_hash, 2, delimiter="#")

    assert params.classes.locations.world.ppl == 4.0


@pytest.mark.unit
def test_override_param(params):
    assert params.classes.sex_types.HM.cis_trans == "cis"
    param_path_hash = "classes#sex_types#HM#cis_trans"

    utils.override_param(params, param_path_hash, "trans", delimiter="#")
    assert params.classes.sex_types.HM.cis_trans == "trans"

    param_path_pipe = "classes|sex_types|HM|cis_trans"
    utils.override_param(params, param_path_pipe, "other")
    assert params.classes.sex_types.HM.cis_trans == "other"

    # test if the override works with an int key
    assert params.partnership.sex.frequency.Sex.bins[1].prob == 0.5
    param_path_haart = "partnership|sex|frequency|Sex|bins|1|prob"
    utils.override_param(params, param_path_haart, 0.0)
    assert params.partnership.sex.frequency.Sex.bins[1].prob == 0.0

    # test that the try/except doesn't silence real key errors
    param_path_fake = "model|time|0"
    with pytest.raises(KeyError):
        utils.override_param(params, param_path_fake, 0)


@pytest.mark.unit
def test_get_independent_bin(params):
    bin_def = params.partnership.sex.frequency.Sex.bins
    rand_gen = FakeRandom(-0.1)

    assert utils.get_independent_bin(rand_gen, bin_def) == 1

    rand_gen = FakeRandom(1.1)
    assert utils.get_independent_bin(rand_gen, bin_def) == len(bin_def)


@pytest.mark.unit
def test_get_cumulative_bin(params):
    bin_def = params.partnership.sex.frequency.Sex.bins
    rand_gen = FakeRandom(-0.1)

    assert utils.get_independent_bin(rand_gen, bin_def) == 1

    rand_gen = FakeRandom(1.1)
    assert utils.get_independent_bin(rand_gen, bin_def) == len(bin_def)


@pytest.mark.unit
def test_grid_file_to_edges():
    grid_file = "tests/params/geography.csv"

    # no diagonals
    edges = utils.grid_file_to_edges(grid_file)

    edge_set = set((e["location_1"], e["location_2"]) for e in edges.values())

    assert ("location_1", "location_3") in edge_set
    assert ("location_1", "location_2") in edge_set
    assert ("location_2", "location_3") in edge_set
    assert ("location_3", "location_4") in edge_set
    assert len(edge_set) == 4

    # with diagonals
    edges = utils.grid_file_to_edges(grid_file, diagonal_neighbors=True)

    edge_set = set((e["location_1"], e["location_2"]) for e in edges.values())

    assert ("location_1", "location_3") in edge_set
    assert ("location_1", "location_2") in edge_set
    assert ("location_2", "location_3") in edge_set
    assert ("location_3", "location_4") in edge_set
    assert ("location_2", "location_4") in edge_set
    assert len(edge_set) == 5


@pytest.mark.unit
def test_grid_file_to_edge_yml(tmpdir):
    grid_file = "tests/params/geography.csv"

    # no diagonals
    outfile_path = os.path.join(tmpdir, "edges.yml")

    utils.grid_file_to_edge_yml(grid_file, outfile_path)

    assert os.path.isfile(outfile_path)

    with open(outfile_path, "r") as f:
        edges = yaml.safe_load(f)

    assert len(edges["edges"]) == 4

    # with diagonals
    outfile_path = os.path.join(tmpdir, "edges_diagonal.yml")

    utils.grid_file_to_edge_yml(grid_file, outfile_path, diagonal_neighbors=True)

    assert os.path.isfile(outfile_path)

    with open(outfile_path, "r") as f:
        edges = yaml.safe_load(f)

    assert len(edges["edges"]) == 5


@pytest.mark.unit
def test_grid2edge_script(tmpdir, monkeypatch):
    grid_file = "tests/params/geography.csv"

    # no diagonals
    outfile_path = os.path.join(tmpdir, "edges.yml")

    monkeypatch.setattr(
        "sys.argv", ["grid2edges", "-i", f"{grid_file}", "-o", f"{outfile_path}"]
    )
    utils.grid_to_edges_cli()

    assert os.path.isfile(outfile_path)

    with open(outfile_path, "r") as f:
        edges = yaml.safe_load(f)

    assert len(edges["edges"]) == 4

    # with diagonals
    outfile_path = os.path.join(tmpdir, "edges_diagonal.yml")

    monkeypatch.setattr(
        "sys.argv", ["grid2edges", "-i", f"{grid_file}", "-o", f"{outfile_path}", "-d"]
    )
    utils.grid_to_edges_cli()

    assert os.path.isfile(outfile_path)

    with open(outfile_path, "r") as f:
        edges = yaml.safe_load(f)

    assert len(edges["edges"]) == 5
