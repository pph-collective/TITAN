import pytest

from conftest import FakeRandom

from titan.location import Location


@pytest.mark.unit
def test_vaccine_init_agent(make_population):
    pop = make_population(n=0)
    pop.params.hiv.start_time = 10
    loc = pop.geography.locations["world"]
    loc.params.vaccine.on_init = True
    pop.pop_random = FakeRandom(-0.1)

    a = pop.create_agent(loc, "white", 0)
    assert a.vaccine.active
    assert a.vaccine.time == 0
    assert a.vaccine.type != ""

    pop.pop_random = FakeRandom(1.0)
    b = pop.create_agent(loc, "white", 0)
    assert b.vaccine.active is False
    assert b.vaccine.time is None
    assert b.vaccine.type == ""


@pytest.mark.unit
def test_vaccine_update_agent(make_model, make_agent):
    model = make_model()
    model.time = model.params.vaccine.start_time
    model.run_random = FakeRandom(-0.1)
    a = make_agent()

    a.prep.active = True
    a.vaccine.update_agent(model)
    assert a.vaccine.active is False

    a.prep.active = False
    a.hiv.active = True
    a.vaccine.update_agent(model)
    assert a.vaccine.active is False

    a.hiv.active = False
    a.vaccine.update_agent(model)
    assert a.vaccine.active is True
    assert a.vaccine.time == model.time

    model.time += (
        a.location.params.demographics[a.race]
        .sex_type[a.sex_type]
        .vaccine.booster.interval
    )
    a.vaccine.update_agent(model)
    assert a.vaccine.active
    assert a.vaccine.time == model.time


@pytest.mark.unit
def test_vaccine_cquisition_risk_multiplier(make_agent):
    a = make_agent()
    assert a.vaccine.get_acquisition_risk_multiplier(0, "sex") == 1.0

    a.vaccine.vaccinate(0)

    assert a.vaccine.get_acquisition_risk_multiplier(1, "sex") < 1.0


@pytest.mark.unit
def test_vaccine_HVTN(make_agent):
    a = make_agent()
    print(a.location.params.vaccine.type)
    a.location.params.vaccine.type = "HVTN702"
    assert a.vaccine.get_acquisition_risk_multiplier(0, "sex") == 1.0

    a.vaccine.vaccinate(0)
    assert a.vaccine.get_acquisition_risk_multiplier(1, "sex") < 1.0


@pytest.mark.unit
def test_vaccine_other(make_agent):
    a = make_agent()
    print(a.location.params.vaccine.type)
    a.location.params.vaccine.type = "other"
    assert a.vaccine.get_acquisition_risk_multiplier(0, "sex") == 1.0

    a.vaccine.vaccinate(0)
    # no risk multiplier
    assert a.vaccine.get_acquisition_risk_multiplier(1, "sex") == 1.0

    a.location.params.vaccine.efficacy = 0.5
    assert a.vaccine.get_acquisition_risk_multiplier(1, "sex") == 0.5
