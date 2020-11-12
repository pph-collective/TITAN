import pytest


@pytest.mark.skip
def test_sex_transmission_do_nothing(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    assert Sex.interact(model, rel) is False

    a.hiv.active = True
    p.hiv.active = True

    # test nothing happens
    assert Sex.interact(model, rel) is False


@pytest.mark.skip
def test_hiv_convert(make_model, make_agent):
    model = make_model()
    a = make_agent()
    a.prep.active = True

    model.run_random = FakeRandom(-0.1)

    model.hiv_convert(a)

    assert a.hiv
    assert a.hiv.time == model.time
    assert a in model.pop.hiv_agents.members
    assert a.prep.active is False


@pytest.mark.skip
def test_diagnose_hiv(make_model, make_agent):
    model = make_model()
    model.params.partner_tracing.prob = 1.0
    model.time = 1
    a = make_agent()
    p = make_agent()
    p.hiv.active = True
    a.partners["Sex"].add(p)

    model.run_random = FakeRandom(1.1)  # always greater than param
    model.diagnose_hiv(a)

    assert a.hiv.dx is False
    assert p.hiv.dx is False
    assert not p.partner_traced

    model.run_random = FakeRandom(-0.1)  # always less than param
    model.diagnose_hiv(a)

    assert a.hiv.dx
    assert a.hiv.dx_time == model.time
    assert p in a.get_partners()
    assert p.partner_traced
    assert p.trace_time == model.time

    assert p.hiv.dx is False
    model.params.demographics[p.race][p.sex_type].hiv.dx.prob = 0

    model.time = p.partner_traced + 1
    model.diagnose_hiv(p)
    assert p.hiv.dx
    assert p.partner_traced is False


@pytest.mark.skip
def test_diagnose_hiv_already_tested(make_model, make_agent):
    model = make_model()
    a = make_agent()

    a.hiv.dx = True

    model.run_random = FakeRandom(-0.1)  # always less than param
    model.diagnose_hiv(a)

    assert a.hiv.dx


@pytest.mark.skip
def test_progress_to_aids_error(make_agent, make_model):
    a = make_agent()
    a.hiv.active = False
    model = make_model()
    a.target_partners = {bond: 0 for bond in model.params.classes.bond_types.keys()}
    model.pop.add_agent(a)
    num_aids = sum([1 for agent in model.pop.hiv_agents if agent.aids])  # get baseline

    # test error case, agent must be HIV+
    with pytest.raises(AssertionError):
        model.progress_to_aids(a)

    assert sum([1 for agent in model.pop.hiv_agents if agent.aids]) == num_aids


@pytest.mark.skip
def test_progress_to_aids_nothing(make_agent, make_model):
    a = make_agent()
    a.hiv.active = True
    model = make_model()
    a.target_partners = {bond: 0 for bond in model.params.classes.bond_types.keys()}
    model.pop.add_agent(a)
    num_aids = sum([1 for agent in model.pop.hiv_agents if agent.aids])  # get baseline

    # test nothing case
    a.hiv.active = True
    a.haart.adherence = 1  # .0051 prob

    model.run_random = FakeRandom(0.9)  # no AIDS

    assert model.progress_to_aids(a) is None
    assert sum([1 for agent in model.pop.hiv_agents if agent.aids]) == num_aids
    assert a.aids is False


@pytest.mark.skip
def test_progress_to_aids_progress(make_agent, make_model):
    a = make_agent()
    a.hiv.active = True
    model = make_model()
    a.target_partners = {bond: 0 for bond in model.params.classes.bond_types.keys()}
    model.pop.add_agent(a)
    num_aids = sum([1 for agent in model.pop.hiv_agents if agent.aids])  # get baseline
    a.location.params.hiv.aids.prob = 1.0

    a.hiv.active = True
    a.haart.adherence = 1  # .0051 prob

    # test progress case
    model.run_random = FakeRandom(0.001)  # AIDS

    assert model.progress_to_aids(a) is None
    assert sum([1 for agent in model.pop.hiv_agents if agent.aids]) == num_aids + 1
    assert a.aids is True


@pytest.mark.skip
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

    assert model.get_transmission_probability("injection", a, p) == p_needle * scale
    assert model.get_transmission_probability("sex", a, p) == p_sex * scale

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

    assert model.get_transmission_probability("sex", a, p) == p_sex_ins * scale
    assert model.get_transmission_probability("sex", p, a) == p_sex_rec * scale


@pytest.mark.skip
def test_get_acute_status(make_agent, params):
    a = make_agent()  # no HIV on init
    assert a.get_acute_status(2) is False
    a.hiv.active = True
    a.hiv.time = 1  # manually force this to test logic
    assert a.get_acute_status(2) is True
