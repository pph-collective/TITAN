import pytest

from titan.interactions import Sex
from titan.agent import Relationship

from conftest import FakeRandom
from titan.parse_params import ObjMap


@pytest.mark.unit
def test_sex_transmission(make_model, make_agent):
    model = make_model()
    model.time = model.params.hiv.start_time
    a = make_agent()
    a.sex_role = "insertive"
    p = make_agent()
    p.sex_role = "receptive"
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    a.hiv = True
    a.hiv_time = model.time  # acute

    rel.total_sex_acts = 0
    model.params.calibration.acquisition = 10

    model.params.calibration.acquisition = 5
    model.params.calibration.sex.act = 10
    model.run_random = FakeRandom(0.6)
    a.location.params.partnership.sex.frequency = ObjMap(
        {"Sex": {"type": "bins", "bins": {1: {"prob": 1.0, "min": 10, "max": 37}}}}
    )
    p.location.params.partnership.sex.frequency = ObjMap(
        {"Sex": {"type": "bins", "bins": {1: {"prob": 1.0, "min": 10, "max": 37}}}}
    )
    # test partner becomes
    Sex.interact(model, rel)
    assert p.hiv

    p.hiv = False

    a.location.params.partnership.sex.frequency = (
        p.location.params.partnership.sex.frequency
    ) = ObjMap(
        {
            "Sex": {
                "type": "distribution",
                "distribution": {
                    "dist_type": "set_value",
                    "vars": {1: {"value": 1, "value_type": "int"}},
                },
            }
        }
    )
    Sex.interact(model, rel)
    assert p.hiv
    # before hiv start time
    p.hiv = False
    model.time = model.params.hiv.start_time - 1
    assert Sex.interact(model, rel) is False
    assert not p.hiv


@pytest.mark.unit
def test_sex_transmission_do_nothing(make_model, make_agent):
    model = make_model()
    model.time = model.params.hiv.start_time
    a = make_agent()
    p = make_agent()
    p_inj = make_agent()
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel_Sex = Relationship(a, p, 10, bond_type="Sex")
    rel_Inj = Relationship(a, p_inj, 10, bond_type="Inj")

    assert Sex.interact(model, rel_Sex) is False

    a.hiv = True
    p.hiv = True

    # test nothing happens
    assert Sex.interact(model, rel_Sex) is False

    a.location.params.partnership.sex.frequency = (
        p.location.params.partnership.sex.frequency
    ) = ObjMap(
        {
            "Sex": {
                "type": "distribution",
                "distribution": {
                    "dist_type": "set_value",
                    "vars": {1: {"value": 0, "value_type": "int"}},
                },
            }
        }
    )
    assert Sex.interact(model, rel_Sex) is False
