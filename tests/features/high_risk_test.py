import pytest

from titan.features import HighRisk

from conftest import FakeRandom


@pytest.mark.unit
def test_become_high_risk(make_model, make_agent):
    model = make_model()
    a = make_agent()

    a.high_risk.become_high_risk(model.pop, model.time, 10)

    assert a.high_risk.active
    assert a.high_risk.ever
    assert a.high_risk.duration == 10
    assert a.high_risk.time == model.time

    a.location.params.features.high_risk = False
    assert not a.high_risk.become_high_risk(model.pop, model.time, 10)


@pytest.mark.unit
def test_update_high_risk(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # try to update when not high risk nor incar
    assert a.high_risk.update_agent(model) is None

    a.high_risk.active = True
    a.high_risk.ever = True
    a.high_risk.duration = 1

    a.high_risk.update_agent(model)

    assert a.high_risk.active
    assert a.high_risk.duration == 0

    a.high_risk.update_agent(model)

    assert a.high_risk.active is False
    assert a.high_risk.ever is True


@pytest.mark.unit
def test_update_high_risk(make_model, make_agent, make_relationship):
    model = make_model()
    model.run_random = FakeRandom(-0.1)
    a = make_agent()
    p = make_agent()
    r = make_relationship(a, p)

    # try to update when not high risk and not incar
    assert a.high_risk.update_agent(model) is None

    # agent incarcerated last time step, partners become high risk
    a.incar.active = True
    a.incar.time = model.time - 1
    a.incar.release_time = model.time

    a.high_risk.update_agent(model)

    assert p.high_risk.active

    # agent released, becomes high risk
    model.time += 1
    a.high_risk.update_agent(model)

    assert a.high_risk.active
    assert a.high_risk.duration == 6
