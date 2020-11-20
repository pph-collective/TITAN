import pytest

from conftest import FakeRandom

from titan.exposures import HIV
from titan.agent import Relationship


@pytest.mark.unit
def test_hiv_expose(make_model, make_agent):
    model = make_model()
    model.run_random = FakeRandom(-0.1)  # forces conversion even if prob is 0
    a = make_agent()
    p = make_agent()
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    HIV.expose(model, "sex", rel, 10)

    assert a.hiv.active is False
    assert p.hiv.active is False

    a.hiv.active = True
    p.hiv.active = True

    # test nothing happens
    HIV.expose(model, "sex", rel, 10)

    assert a.hiv.active
    assert p.hiv.active

    p.hiv.active = False

    # test conversion happens
    HIV.expose(model, "sex", rel, 10)

    assert a.hiv.active
    assert p.hiv.active


@pytest.mark.unit
def test_hiv_init(make_population, make_agent):
    pop = make_population()
    pop.pop_random = FakeRandom(-0.1)
    a = make_agent()

    time = pop.params.hiv.start_time - 1

    a.hiv.init_agent(pop, time)

    assert a.hiv.active is False

    time = pop.params.hiv.start_time

    a.hiv.init_agent(pop, time)

    assert a.hiv.active
    assert a.hiv.time == time - pop.params.hiv.max_init_time
    assert a.hiv.dx
    assert a.hiv.dx_time == a.hiv.time
    assert a.hiv.aids
    assert a in HIV.agents
    assert HIV.dx_counts[a.race][a.sex_type] == 1


@pytest.mark.unit
def test_hiv_convert(make_model, make_agent):
    model = make_model()
    a = make_agent()
    a.prep.active = True

    model.run_random = FakeRandom(-0.1)

    a.hiv.convert(model)

    assert a.hiv.active
    assert a.hiv.time == model.time
    assert a in HIV.agents
    assert a.prep.active is False


@pytest.mark.unit
def test_diagnose_hiv(make_model, make_agent):
    model = make_model()
    model.params.partner_tracing.prob = 1.0
    model.time = 1
    a = make_agent()
    p = make_agent()
    p.hiv.active = True
    a.partners["Sex"].add(p)

    model.run_random = FakeRandom(-0.1)  # always less than param
    a.hiv.diagnose(model)

    assert a.hiv.dx
    assert a.hiv.dx_time == model.time


@pytest.mark.unit
def test_get_transmission_probability(make_model, make_agent):
    model = make_model()
    a = make_agent(race="white", SO="MSM")
    a.haart.adherence = 1
    a.sex_role = "versatile"

    p = make_agent(race="white", SO="MSM")
    p.sex_role = "versatile"
    p.haart.adherence = 1

    # test versatile-versatile relationship
    p_needle = (
        model.params.partnership.injection.transmission.haart_scaling[1].scale
        * model.params.partnership.injection.transmission.base
    )
    p_sex = (
        model.params.partnership.sex.haart_scaling["MSM"][1].prob
        * model.params.partnership.sex.acquisition.MSM.versatile
    )
    scale = model.params.calibration.acquisition

    assert (
        a.hiv.get_transmission_probability(model, "injection", p, 1) == p_needle * scale
    )
    assert a.hiv.get_transmission_probability(model, "sex", p, 1) == p_sex * scale

    # test one vers agent, one receptive agent
    a.sex_role = "receptive"
    p_sex_ins = (
        model.params.partnership.sex.haart_scaling.MSM[1].prob
        * model.params.partnership.sex.acquisition.MSM.insertive
    )
    p_sex_rec = (
        model.params.partnership.sex.haart_scaling.MSM[1].prob
        * model.params.partnership.sex.acquisition.MSM.receptive
    )

    assert a.hiv.get_transmission_probability(model, "sex", p, 1) == p_sex_ins * scale
    assert p.hiv.get_transmission_probability(model, "sex", a, 1) == p_sex_rec * scale


@pytest.mark.unit
def test_get_acute_status(make_agent, make_model):
    model = make_model()
    a = make_agent()  # no HIV on init
    assert a.hiv.get_acute_status(model.time + 2) is False
    a.hiv.convert(model)
    assert a.hiv.get_acute_status(model.time + 2) is True
