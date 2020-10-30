import pytest

from titan.interactions import Sex
from titan.agent import Relationship

from conftest import FakeRandom


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
    a.hiv.time = model.time  # acute

    rel.total_sex_acts = 0
    model.params.calibration.acquisition = 10

    model.params.calibration.acquisition = 5
    model.params.calibration.sex.act = 10
    model.run_random = FakeRandom(0.6)

    # test partner becomes
    Sex.interact(model, rel)
    assert p.hiv


@pytest.mark.unit
def test_sex_transmission_do_nothing(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    assert Sex.interact(model, rel) is False

    a.hiv = True
    p.hiv = True

    # test nothing happens
    assert Sex.interact(model, rel) is False
