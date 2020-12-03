import pytest

from conftest import FakeRandom


@pytest.mark.unit
def test_update_haart_t1(make_model, make_agent):
    model = make_model()
    model.time = 1
    a = make_agent(race="white")

    a.hiv = True

    # nothing happens, not tested
    a.haart.update_agent(model)
    assert a.haart.adherent is False
    assert a.haart.active is False

    # t0 agent initialized HAART
    a.hiv_dx = True

    # go on haart
    model.run_random = FakeRandom(
        -0.1
    )  # means this will always be less than params even though not possible in reality
    a.haart.update_agent(model)

    assert a.haart.adherent is True
    assert a.haart.active

    # go off haart
    a.haart.update_agent(model)

    assert a.haart.adherent is False
    assert a.haart.active is False

    # Try haart cap with low cap, one agent on haart and one diagnosed agent. Nothing happens
    a.location.params.hiv.haart_cap = True
    a.haart.counts[a.race][a.sex_type] = 1
    a.location.params.demographics[a.race][a.sex_type].haart.prob = 0.1
    model.pop.dx_counts[a.race][a.sex_type] = 1
    a.haart.update_agent(model)
    assert a.haart.active is False
    assert a.haart.adherent is False

    # Increase cap. Agent goes on prep
    a.location.params.demographics[a.race][a.sex_type].haart.prob = 2.0
    a.haart.update_agent(model)
    assert a.haart.active
    assert a.haart.adherent is True
