import pytest
import os
import shutil
import csv
import numpy as np

import titan.simulation_lib as sl


@pytest.fixture
def setup_results_dir():
    outfile_dir = os.path.join(os.getcwd(), "results")
    if os.path.isdir(outfile_dir):
        shutil.rmtree(outfile_dir)
    os.mkdir(outfile_dir)
    yield outfile_dir
    shutil.rmtree(outfile_dir)


def test_initiate_result_dict():
    # initiate result dict with 2 time steps
    d = sl.initiate_result_dict(2)

    assert "Prv_HIV" in d
    assert 2 in d["Prv_HIV"]
    assert d["Prv_HIV"][2] == []


def test_safe_divide():
    assert sl.safe_divide(0, 0) == 0.0
    assert sl.safe_divide(0, 1) == 0.0
    assert sl.safe_divide(1, 2) == 0.5


@pytest.mark.integration
def test_simulation_one_rep(setup_results_dir):
    result_dict = sl.simulation(1, 2, 10, 2, 2, 2)

    assert "Prv_HIV" in result_dict
    assert 2 in result_dict["Prv_HIV"]
    assert len(result_dict["Prv_HIV"][2]) == 1


@pytest.mark.integration
def test_simulation_multi_rep(setup_results_dir):
    result_dict = sl.simulation(3, 2, 10, 2, 2, 2)

    assert "Prv_HIV" in result_dict
    assert 2 in result_dict["Prv_HIV"]
    assert len(result_dict["Prv_HIV"][2]) == 3


@pytest.mark.integration
def test_save_results(tmpdir, setup_results_dir):
    result_dict = sl.simulation(3, 2, 10, 2, 2, 2)
    dir_str = str(tmpdir)
    sim_num = 123
    sl.save_results(2, result_dict, dir_str, sim_num)

    # check the file exists
    result_file = os.path.join(dir_str, f"Result_simulation_{sim_num}.txt")
    assert os.path.isfile(result_file)

    # spot check the values match as expected
    with open(result_file, newline="") as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            t = i + 1
            print(t)
            print(result_dict["n_Relations"][t])
            print(row["n_Relations_mean"])
            num_rels_mean = np.mean(result_dict["n_Relations"][t])
            assert round(num_rels_mean, 5) == float(row["n_Relations_mean"])

