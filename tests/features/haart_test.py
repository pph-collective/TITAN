import pytest

from conftest import FakeRandom


@pytest.mark.unit
def test_update_haart_t1(make_model, make_agent):
    model = make_model()
    model.time = 1
    a = make_agent(race="white")

    a.hiv.active = True

    # nothing happens, not tested
    a.haart.update_agent(model)
    assert a.haart.adherent is False
    assert a.haart.active is False

    # t0 agent initialized HAART
    a.hiv.dx = True
    a.hiv.dx_time = model.time

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
    a.hiv.dx_counts[a.race][a.sex_type] = 5
    a.location.params.demographics[a.race].sex_type[a.sex_type].drug_type[
        a.drug_type
    ].haart.prob = 0.1
    a.haart.update_agent(model)
    assert a.haart.active is False
    assert a.haart.adherent is False

    # Increase cap. Agent goes on haart
    a.location.params.demographics[a.race].sex_type[a.sex_type].drug_type[
        a.drug_type
    ].haart.cap = 2.0
    a.haart.update_agent(model)
    assert a.haart.active
    assert a.haart.adherent is True

    # falls off adherence
    model.run_random = FakeRandom(1.0)
    a.location.params.demographics[a.race].sex_type[a.sex_type].drug_type[
        a.drug_type
    ].haart.adherence.discontinue = 5.0
    a.haart.update_agent(model)
    assert a.haart.active
    assert a.haart.adherent is False


def test_update_haart_dx_time(make_model, make_agent):
    model = make_model()
    model.time = 1
    a = make_agent(race="white")
    model.run_random = FakeRandom(0.5)

    a.hiv.active = True
    a.hiv.dx = True
    a.hiv.dx_time = 1  # agent dx duration is 0
    a.haart.update_agent(model)
    assert not a.haart.active

    a.hiv.dx_time -= 10
    a.haart.update_agent(model)
    assert a.haart.active


def test_update_haart_reinit(make_model, make_agent):
    model = make_model()
    a = make_agent(race="white")

    a.hiv.dx = True
    a.hiv.dx_time = model.time - 1  # make duration 1, would pass if not reinit
    a.haart.ever = True  # uses reinit prob
    a.haart.update_agent(model)
    assert not a.haart.active
