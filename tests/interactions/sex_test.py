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

    a.hiv.active = True
    a.hiv.time = model.time  # acute

    rel.total_sex_acts = 0
    model.params.calibration.acquisition = 10

    model.params.calibration.acquisition = 5
    model.params.calibration.sex.act = 10
    model.run_random = FakeRandom(0.6)

    # test partner becomes
    Sex.interact(model, rel)
    assert p.hiv.active


@pytest.mark.unit
def test_sex_num_acts(make_model, make_agent, make_relationship, params):
    params.hiv.dx.risk_reduction.sex = 1.0
    model = make_model()
    model.np_random = FakeRandom(1.0)
    a = make_agent()
    p = make_agent()
    rel = make_relationship(a, p)

    assert Sex.get_num_acts(model, rel) > 0

    a.hiv.active = True
    a.hiv.dx = True

    assert Sex.get_num_acts(model, rel) == 0
