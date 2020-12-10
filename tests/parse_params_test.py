import pytest
import os

import titan.parse_params as pp


setting = "tests/params/setting.yml"


@pytest.mark.unit
def test_create_params_setting(tmpdir):
    param_file_clean = "tests/params/setting_params.yml"
    params = pp.create_params(setting, param_file_clean, tmpdir)

    result_file = os.path.join(tmpdir, "params.yml")
    assert os.path.isfile(result_file)

    assert "HM" not in list(params.classes.sex_types.keys())
    assert "white" not in params.classes.races
    assert params.partnership.sex.acquisition.MSM.versatile == 0.00745  # from setting
    assert (
        params.partnership.sex.acquisition.MSM.insertive == 0.0011
    )  # from param, rounded


@pytest.mark.unit
def test_create_params_error(tmpdir):
    param_file_error = "tests/params/basic.yml"
    # should error since base has demographics for more populations than there are classes in settings, so the ppl's don't add up
    with pytest.raises(AssertionError):
        pp.create_params(setting, param_file_error, tmpdir)

    # this params file has a bad value, this makes sure the correct key is included in the error message
    param_file_bad_val = "tests/params/basic_seeded_error.yml"
    with pytest.raises(AssertionError) as excinfo:
        pp.create_params(setting, param_file_bad_val, tmpdir)

    assert (
        "[.demographics.black.sex_type.MSM.drug_type.Inj.num_partners.Sex.dist_type]"
        in str(excinfo.value)
    )
