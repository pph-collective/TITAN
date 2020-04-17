import pytest

import os
import subprocess
import csv


@pytest.mark.integration
def test_model_runs():
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic.yml"
    )

    subprocess.check_call([f, f"-p {param_file}"])
    assert True


@pytest.mark.integration
def test_model_reproducible(tmpdir):
    path_a = tmpdir.mkdir("result_a")
    path_b = tmpdir.mkdir("result_b")
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "basic_seeded.yml"
    )

    subprocess.check_call([f, f"-p {param_file}", f"-o {path_a}"])
    subprocess.check_call([f, f"-p {param_file}", f"-o {path_b}"])

    result_file_a = os.path.join(path_a, "basicReport_MSM_BLACK.txt")
    result_file_b = os.path.join(path_b, "basicReport_MSM_BLACK.txt")
    assert os.path.isfile(result_file_a)
    with open(result_file_a, newline="") as fa, open(result_file_b, newline="") as fb:
        reader_a = csv.DictReader(fa, delimiter="\t")
        res_a = []
        for row in reader_a:
            res_a.append(row)

        reader_b = csv.DictReader(fb, delimiter="\t")
        res_b = []
        for row in reader_b:
            res_b.append(row)

    for i in range(len(res_a)):
        assert res_a[i]["t"] == res_b[i]["t"]
        assert res_a[i]["run_id"] != res_b[i]["run_id"]
        assert res_a[i]["rseed"] == res_b[i]["rseed"]
        assert res_a[i]["pseed"] == res_b[i]["pseed"]
        assert res_a[i]["Total"] == res_b[i]["Total"]
        assert res_a[i]["HIV"] == res_b[i]["HIV"]
        assert res_a[i]["PrEP"] == res_b[i]["PrEP"]
        assert res_a[i]["Deaths"] == res_b[i]["Deaths"]
        assert res_a[i]["HIV"] == res_b[i]["HIV"]


@pytest.mark.integration
def test_model_settings_run(tmpdir):
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "run_titan.py")
    param_file = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "params", "integration_base.yml"
    )

    for item in os.listdir(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "settings")
    ):
        if "__" not in item and item != "base":
            path = tmpdir.mkdir(item)
            subprocess.check_call([f, f"-p {param_file}", f"-o {path}", f"-S {item}"])
            assert True
