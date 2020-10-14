import pytest

from titan.features import HighRisk


@pytest.mark.unit
def test_become_high_risk(make_model, make_agent):
    a = make_agent()

    a.high_risk.become_high_risk(1, 10)

    assert a.high_risk.active
    assert a.high_risk.ever
    assert a.high_risk.duration == 10
    assert a.high_risk.time == 1

    a.location.params.features.high_risk = False
    assert not a.high_risk.become_high_risk(1, 10)


@pytest.mark.unit
def test_update_high_risk(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # try to update when not high risk
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
