import pytest
import os

import titan.parse_params as pp


setting = "tests/params/setting.yml"


def test_create_params_base(tmpdir):
    param_file_clean = "tests/params/setting_params.yml"
    params = pp.create_params(setting, param_file_clean, tmpdir, use_base=True)

    result_file = os.path.join(tmpdir, "params.yml")
    assert os.path.isfile(result_file)

    assert "HM" not in list(params.classes.sex_types.keys())
    assert "WHITE" not in params.classes.races
    assert params.partnership.sex.acquisition.MSM.versatile == 0.00745  # from setting
    assert params.partnership.sex.acquisition.MSM.receptive == 0.0138  # from base
    assert params.partnership.sex.acquisition.MSM.insertive == 0.0011  # from param
    assert params.demographics.BLACK.MSM.num_partners.rand == 1.0


def test_create_params_error(tmpdir):
    param_file_error = "tests/params/basic.yml"
    # should error since base has demographics for more populations than there are classes in settings, so the ppl's don't add up
    with pytest.raises(AssertionError):
        pp.create_params(setting, param_file_error, tmpdir)


def test_check_item_min_max():
    defs = {"min": 0, "max": 3, "type": "float"}

    # happy case
    assert pp.check_item(1.5, defs) == 1.5

    # type coerciion to float
    assert pp.check_item(1, defs) == 1.0

    # below min
    with pytest.raises(AssertionError):
        pp.check_item(-1.5, defs)

    # above max
    with pytest.raises(AssertionError):
        pp.check_item(4.5, defs)


def test_check_item_int():
    defs = {"min": 0, "max": 3, "type": "int"}

    # happy case
    assert pp.check_item(1, defs) == 1

    # not int
    with pytest.raises(AssertionError):
        pp.check_item(1.5, defs)


def test_check_item_boolean():
    defs = {"type": "boolean"}

    # happy case
    assert pp.check_item(True, defs) == True

    # not bool
    with pytest.raises(AssertionError):
        pp.check_item(1, defs)


def test_check_item_enum():
    defs = {"type": "enum", "values": ["a", "b"]}

    # happy case
    assert pp.check_item("a", defs) == "a"

    # not bool
    with pytest.raises(AssertionError):
        pp.check_item("c", defs)


def test_check_item_array():
    defs = {"type": "enum", "values": ["a", "b"]}

    # happy case
    assert pp.check_item(["a", "b"], defs) == ["a", "b"]
    assert pp.check_item([], defs) == []
    assert pp.check_item(["b"], defs) == ["b"]

    # not in array
    with pytest.raises(AssertionError):
        pp.check_item(["c"], defs)

    # some not in array
    with pytest.raises(AssertionError):
        pp.check_item(["a", "c"], defs)


def test_check_item_array():
    defs = {"type": "keys"}

    # happy case
    assert pp.check_item(["a", "b"], defs, keys=["a", "b"]) == ["a", "b"]
    assert pp.check_item([], defs, keys=["a", "b"]) == []
    assert pp.check_item(["b"], defs, keys=["a", "b"]) == ["b"]

    # not in array
    with pytest.raises(AssertionError):
        pp.check_item(["c"], defs, keys=["a", "b"])

    # some not in array
    with pytest.raises(AssertionError):
        pp.check_item(["a", "c"], defs, keys=["a", "b"])
