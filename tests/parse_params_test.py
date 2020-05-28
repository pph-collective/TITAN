import pytest
import os

import titan.parse_params as pp


setting = "tests/params/setting.yml"


@pytest.mark.unit
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


@pytest.mark.unit
def test_create_params_error(tmpdir):
    param_file_error = "tests/params/basic.yml"
    # should error since base has demographics for more populations than there are classes in settings, so the ppl's don't add up
    with pytest.raises(AssertionError):
        pp.create_params(setting, param_file_error, tmpdir)

    param_file_bad_val = "tests/params/basic_seeded_error.yml"
    with pytest.raises(AssertionError) as excinfo:
        pp.create_params(setting, param_file_bad_val, tmpdir)

        assert "[.demographics.BLACK.MSM.num_partners.Sex.disk_type]" in str(
            excinfo.value
        )


@pytest.mark.unit
def test_check_item_min_max():
    defs = {"min": 0, "max": 3, "type": "float"}
    key_path = "item.test"

    # happy case
    assert pp.check_item(1.5, defs, key_path) == 1.5

    # type coerciion to float
    assert pp.check_item(1, defs, key_path) == 1.0

    # below min
    with pytest.raises(AssertionError):
        pp.check_item(-1.5, defs, key_path)

    # above max
    with pytest.raises(AssertionError):
        pp.check_item(4.5, defs, key_path)


@pytest.mark.unit
def test_check_item_int():
    defs = {"min": 0, "max": 3, "type": "int"}
    key_path = "item.test"

    # happy case
    assert pp.check_item(1, defs, key_path) == 1

    # not int
    with pytest.raises(AssertionError):
        pp.check_item(1.5, defs, key_path)


@pytest.mark.unit
def test_check_item_boolean():
    defs = {"type": "boolean"}
    key_path = "item.test"

    # happy case
    assert pp.check_item(True, defs, key_path) == True

    # not bool
    with pytest.raises(AssertionError):
        pp.check_item(1, defs, key_path)


@pytest.mark.unit
def test_check_item_enum():
    defs = {"type": "enum", "values": ["a", "b"]}
    key_path = "item.test"

    # happy case
    assert pp.check_item("a", defs, key_path) == "a"

    # not bool
    with pytest.raises(AssertionError):
        pp.check_item("c", defs, key_path)


@pytest.mark.unit
def test_check_item_array():
    defs = {"type": "array", "values": ["a", "b"]}
    key_path = "item.test"

    # happy case
    assert pp.check_item(["a", "b"], defs, key_path) == ["a", "b"]
    assert pp.check_item([], defs, key_path) == []
    assert pp.check_item(["b"], defs, key_path) == ["b"]

    # not in array
    with pytest.raises(AssertionError):
        pp.check_item(["c"], defs, key_path)

    # some not in array
    with pytest.raises(AssertionError):
        pp.check_item(["a", "c"], defs, key_path)


@pytest.mark.unit
def test_check_item_keys():
    defs = {"type": "keys"}
    key_path = "item.test"

    # happy case
    assert pp.check_item(["a", "b"], defs, key_path, keys=["a", "b"]) == ["a", "b"]
    assert pp.check_item([], defs, key_path, keys=["a", "b"]) == []
    assert pp.check_item(["b"], defs, key_path, keys=["a", "b"]) == ["b"]

    # not in array
    with pytest.raises(AssertionError):
        pp.check_item(["c"], defs, key_path, keys=["a", "b"])

    # some not in array
    with pytest.raises(AssertionError):
        pp.check_item(["a", "c"], defs, key_path, keys=["a", "b"])


@pytest.mark.unit
def test_check_item_class_enum():
    defs = {"type": "enum", "class": "bond_types"}
    key_path = "item.test"

    nested_pops = {"bond_types": {"Inj": {"a": 1}, "Other": {"b": 2}}}
    flat_pops = {"bond_types": ["Inj", "Other"]}

    # happy case
    assert pp.check_item("Inj", defs, key_path, pops=nested_pops) == "Inj"
    assert pp.check_item("Inj", defs, key_path, pops=flat_pops) == "Inj"

    # not in array
    with pytest.raises(AssertionError):
        pp.check_item("Junk", defs, key_path, pops=nested_pops)

    # some not in array
    with pytest.raises(AssertionError):
        pp.check_item("Junk", defs, key_path, pops=flat_pops)


@pytest.mark.unit
def test_check_item_class_array():
    defs = {"type": "array", "class": "bond_types"}
    key_path = "item.test"

    nested_pops = {"bond_types": {"Inj": {"a": 1}, "Other": {"b": 2}}}
    flat_pops = {"bond_types": ["Inj", "Other"]}

    # happy case
    assert pp.check_item(["Inj"], defs, key_path, pops=nested_pops) == ["Inj"]
    assert pp.check_item(["Inj"], defs, key_path, pops=flat_pops) == ["Inj"]

    # not in array
    with pytest.raises(AssertionError):
        pp.check_item(["Junk"], defs, key_path, pops=nested_pops)

    # some not in array
    with pytest.raises(AssertionError):
        pp.check_item(["Junk"], defs, key_path, pops=flat_pops)
