import pytest

from conftest import FakeRandom

from titan.exposures import MonkeyPox
from titan.agent import Relationship


@pytest.mark.unit
def test_monkeypox_expose(make_model, make_agent):
    model = make_model()
    model.run_random = FakeRandom(0.0)  # always less than param
    a = make_agent()
    p = make_agent()
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    MonkeyPox.expose(model, "sex", rel, 1)

    assert a.monkeypox.active is False
    assert p.monkeypox.active is False

    a.monkeypox.active = True
    a.monkeypox.time = model.time
    p.monkeypox.active = True
    a.monkeypox.time = model.time

    # test nothing happens
    MonkeyPox.expose(model, "sex", rel, 10)

    assert a.monkeypox.active
    assert p.monkeypox.active

    p.monkeypox.active = False

    # test conversion happens
    MonkeyPox.expose(model, "sex", rel, 10)

    assert a.monkeypox.active
    assert p.monkeypox.active


@pytest.mark.unit
def test_monkeypox_convert(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.run_random = FakeRandom(-0.1)

    a.monkeypox.convert(model)

    assert a.monkeypox.active
    assert a.monkeypox.time == model.time
    assert a in MonkeyPox.agents


# TODO talk with Ellen about why this is failing?
# @pytest.mark.unit
# def test_diagnose_monkeypox(make_model, make_agent):
#     model = make_model()
#     model.params.partner_tracing.prob = 1.0
#     model.time = 1
#     a = make_agent()
#     p = make_agent()
#     p.monkeypox.active = True
#     a.partners["Sex"].add(p)

#     model.run_random = FakeRandom(-0.1)  # always less than param
#     a.monkeypox.diagnose(model)

#     assert a.monkeypox.dx
#     assert a.monkeypox.dx_time == model.time


@pytest.mark.unit
def test_get_transmission_probability(make_model, make_agent):
    model = make_model()
    a = make_agent(race="white", SO="MSM")
    a.sex_role = "versatile"

    p = make_agent(race="white", SO="MSM")
    p.sex_role = "versatile"

    p_sex = model.params.partnership.sex.acquisition.MSM.versatile
    scale = model.params.calibration.acquisition

    # must be in acute phase
    a.monkeypox.active = True
    a.monkeypox.time = model.time

    assert a.monkeypox.get_transmission_probability(model, "sex", p, 1) == p_sex * scale


@pytest.mark.unit
def test_get_acute_status(make_agent, make_model):
    model = make_model()
    a = make_agent()
    assert a.monkeypox.get_acute_status(model.time + 2) is False
    a.monkeypox.convert(model)
    assert a.monkeypox.get_acute_status(model.time + 2) is True
