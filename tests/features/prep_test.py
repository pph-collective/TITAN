import pytest

from conftest import FakeRandom

from titan.features import Prep
from titan.agent import Relationship


@pytest.mark.unit
def test_init_agent(make_agent, make_population, make_model):
    model = make_model()
    pop = make_population(n=1)

    # test allcomers (no prep)
    pop.pop_random = FakeRandom(0.9)
    a1 = make_agent(race="black")
    a2 = make_agent(race="white")
    assert not a1.prep.active and not a2.prep.active

    # test allcomers
    pop.pop_random = FakeRandom(0.4)
    a1.prep.init_agent(pop, 10)
    a2.prep.init_agent(pop, 10)
    assert a1.prep.active
    assert a2.prep.active

    # reset prep, change to racial model between black and white init prob
    a1.prep.active = a2.prep.active = False
    a1.location.params.prep.target_model = a2.location.params.prep.target_model = [
        "Racial"
    ]
    a1.prep.init_agent(pop, 10)
    a2.prep.init_agent(pop, 10)
    assert a1.prep.active
    assert not a2.prep.active


@pytest.mark.unit
def test_num_partners(make_agent, make_population, make_model, make_relationship):
    model = make_model()
    model.params.prep.target_model = ["Allcomers", "top_partners"]
    model.params.prep.top_partners = 0.01
    model.pop.pop_random = FakeRandom(0.1)
    # add new agent to population
    a1 = make_agent()
    model.pop.add_agent(a1)
    assert not a1.prep.active

    # make sure this agent has most partners
    for ag in model.pop.all_agents.members:
        if ag != a1:
            make_relationship(ag, a1)
    for feature in model.features:
        feature.update_pop(model)
    a1.prep.init_agent(model.pop, 0)
    assert a1.prep.active
    assert not model.pop.all_agents.members.pop().prep.active


@pytest.mark.unit
def test_discontinue_prep(make_agent):
    a = make_agent()

    # set up so the agent appears to be on PrEP
    a.prep.active = True
    a.prep.time = 0
    a.prep.last_dose_time = 0
    num_prep = Prep.counts[a.race]
    Prep.counts[a.race] += 1

    a.prep.discontinue()

    assert a.prep.active is False
    assert a.prep.time is None
    assert a.prep.last_dose_time is None
    assert num_prep == Prep.counts[a.race]


@pytest.mark.unit
def test_update_agent_prep_decrement_time(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.run_random = FakeRandom(1.1)
    model.time = model.params.prep.start_time

    # set up so the agent appears to be on PrEP
    a.prep.active = True
    a.prep.type = "Oral"
    a.prep.last_dose_time = model.time - 1
    num_prep = Prep.counts[a.race]
    Prep.counts[a.race] += 1

    a.prep.update_agent(model)

    assert a.prep.active
    assert a.prep.last_dose_time == model.time
    assert num_prep + 1 == Prep.counts[a.race]


@pytest.mark.unit
def test_update_agent_prep_decrement_end(make_model, make_agent):
    model = make_model()
    a = make_agent(race="white")

    model.run_random = FakeRandom(-0.1)
    model.time = model.params.prep.start_time

    # set up so the agent appears to be on PrEP
    a.prep.active = True
    a.prep.type = "Oral"
    a.prep.last_dose_time = model.time - 1
    num_prep = Prep.counts[a.race]
    Prep.counts[a.race] += 1

    a.prep.update_agent(model)

    assert a.prep.active is False
    assert a.prep.type == ""
    assert a.prep.last_dose_time is None
    assert num_prep == Prep.counts[a.race]


@pytest.mark.unit
def test_update_agent_prep_decrement_not_end(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.run_random = FakeRandom(1.1)
    model.time = model.params.prep.start_time

    # set up so the agent appears to be on PrEP
    a.prep.active = True
    a.prep.last_dose_time = model.time
    a.prep.type = "Inj"
    num_prep = Prep.counts[a.race]
    Prep.counts[a.race] += 1

    a.prep.update_agent(model)

    assert a.prep.active
    assert a.prep.last_dose_time == model.time
    assert num_prep + 1 == Prep.counts[a.race]


@pytest.mark.unit
def test_update_agent_prep_decrement_inj_end(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.run_random = FakeRandom(1.1)
    model.time = model.params.prep.start_time

    # set up so the agent appears to be on PrEP
    a.prep.active = True
    a.prep.last_dose_time = (
        model.time - 12
    )  # last step before hitting year mark and discontinuing
    a.prep.type = "Inj"
    num_prep = Prep.counts[a.race]
    Prep.counts[a.race] += 1

    a.prep.update_agent(model)

    assert a.prep.active is False
    assert (
        a.prep.last_dose_time is None
    )  # 3 -> -1 -> +1 == 0 # Inj no longer in PrEP types
    assert num_prep == Prep.counts[a.race]


@pytest.mark.unit
def test_initiate_prep_assertions(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # no PreP if already PreP
    a.prep.active = True
    assert a.prep.initiate(model) is None

    # no PrEP if already HIV
    a.prep.active = False
    a.hiv.active = True
    assert a.prep.initiate(model) is None


@pytest.mark.unit
def test_initiate_prep_force_adh(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # forcing, adherant, inj
    model.run_random = FakeRandom(-0.1)
    a.prep.initiate(model, True)
    assert a.prep.active
    assert a.prep.adherent is True
    assert a.prep.last_dose_time == model.time


@pytest.mark.unit
def test_initiate_prep_force_non_adh(make_model, make_agent):
    model = make_model()
    a = make_agent()
    # forcing, non-adherant, inj
    model.run_random = FakeRandom(1.0)
    a.prep.initiate(model, True)
    assert a.prep.active
    assert a.prep.adherent is False
    assert a.prep.last_dose_time == model.time


@pytest.mark.unit
def test_initiate_prep_eligible(make_model, make_agent):
    model = make_model()

    # make sure there's room to add more prep agents
    a = make_agent(SO="HF")  # model is "CDCwomen"
    p = make_agent(DU="Inj")
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    p.hiv.dx = True
    p.external_exposure.active = True
    model.time = 10
    a.location.params.prep.cap = 1.0
    a.location.params.prep.target_model = ["cdc_women"]
    rel = Relationship(a, p, 10, bond_type="Sex")
    # non-forcing, adherant, inj
    model.run_random = FakeRandom(-0.1)
    a.prep.initiate(model)
    assert a.prep.active
    assert a.prep.adherent is True
    assert a.prep.last_dose_time == 10
    assert a.prep.time == 10


@pytest.mark.unit
def test_initiate_prep_eligible_racial(make_model, make_agent):
    model = make_model()

    # make sure there's room to add more prep agents
    a = make_agent(SO="HF")  # model is "CDCwomen"
    p = make_agent(DU="Inj")
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    p.hiv.dx = True
    p.external_exposure.active = True
    model.time = 10
    a.location.params.prep.cap = 1.0
    a.location.params.prep.target_model = ["Racial"]
    rel = Relationship(a, p, 10, bond_type="Sex")
    # non-forcing, adherant, inj
    model.run_random = FakeRandom(-0.1)
    a.prep.initiate(model)
    assert a.prep.active
    assert a.prep.adherent is True
    assert a.prep.last_dose_time == 10
    assert a.prep.time == 10


@pytest.mark.unit
def test_cdc_eligible(make_agent, make_relationship):
    # test MSM
    a = make_agent()
    p = make_agent()
    r = make_relationship(a, p)
    assert a.prep.cdc_eligible()

    # test WSW fail
    a = make_agent(SO="WSW")
    assert not a.is_msm()
    assert not a.prep.cdc_eligible()

    # relationship not eligible
    p = make_agent(SO="HM")
    r = make_relationship(a, p)
    assert not a.prep.cdc_eligible()

    # relationship eligible
    p.hiv.dx = True
    assert a.prep.cdc_eligible()

    # ongoing duration fail
    a.location.params.partnership.ongoing_duration = 10
    assert not a.prep.cdc_eligible()


@pytest.mark.unit
def test_prep_eligible(make_agent, make_relationship):
    a = make_agent(SO="HF")
    p = make_agent(SO="HM")

    # test no model
    a.location.params.prep.target_model = ["cdc_women"]
    assert not a.prep.eligible(10)

    # test cdc_women
    a.location.params.prep.target_model.append("cdc_women")
    assert not a.prep.eligible(10)
    r = make_relationship(a, p)
    assert not p.is_msm()
    assert not a.prep.eligible(10)
    p.drug_type = "Inj"
    assert a.prep.eligible(10)

    # test Allcomers and Racial
    a.location.params.prep.target_model.append("Allcomers")
    assert a.prep.eligible(10)
    a.location.params.prep.target_model = ["Racial"]
    assert a.prep.eligible(10)

    # test cdc_msm
    a.location.params.prep.target_model = ["cdc_msm"]
    assert not a.prep.eligible(10)
    msm_agent = make_agent()
    msm_agent.location.params.prep.target_model = ["cdc_msm"]
    make_relationship(msm_agent, p)
    assert msm_agent.is_msm()
    assert msm_agent.prep.eligible(10)

    # test pwid sex
    a.location.params.prep.target_model = ["pwid_sex"]
    p.location.params.prep.target_model = ["pwid_sex"]
    assert not a.prep.eligible(10)
    assert p.prep.eligible(10)

    # test pwid
    a.location.params.prep.target_model = ["pwid"]
    p.location.params.prep.target_model = ["pwid"]
    assert not a.prep.eligible(10)
    assert p.prep.eligible(10)
    p = make_agent(DU="Inj", SO="HM")
    p.location.params.prep.target_model = ["pwid"]
    assert p.prep.eligible(10)  # still eligible without partners

    # test ssp
    a.location.params.prep.target_model = ["ssp"]
    p.location.params.prep.target_model = ["ssp"]
    assert not a.prep.eligible(10)
    assert not p.prep.eligible(10)
    p.syringe_services.active = True
    assert p.prep.eligible(10)

    p.location.params.prep.target_model = ["ssp_sex"]
    assert not p.prep.eligible(10)
    make_relationship(p, msm_agent)
    assert p.prep.eligible(10)

    b = make_agent()
    b.location.params.prep.target_model = ["Allcomers"]
    assert b.prep.eligible(10)
    b.prep.active = True
    assert not b.prep.eligible(10)
    b.prep.active = False
    b.hiv.active = True
    assert not b.prep.eligible(10)
    b.hiv.active = False
    b.location.params.features.random_trial = True
    assert not b.prep.eligible(10)


@pytest.mark.unit
def test_enroll_prep_choice(make_agent, params):
    rand_gen = FakeRandom(-0.1)
    a = make_agent()
    a.location.params.prep.type = ["Oral", "Inj"]

    a.prep.enroll(rand_gen, 0)

    assert a.prep.active
    assert a.prep.last_dose_time == 0
    assert a.prep.adherent is True
    assert a.prep.type == "Inj"


@pytest.mark.unit
def test_enroll_prep_one(make_agent, params):
    rand_gen = FakeRandom(1.1)
    a = make_agent()
    a.location.params.prep.type = ["Oral"]

    a.prep.enroll(rand_gen, 0)

    assert a.prep.active
    assert a.prep.last_dose_time == 0
    assert a.prep.adherent is False
    assert a.prep.type == "Oral"


@pytest.mark.unit
def test_progress_inj_prep(make_agent, params, make_model):
    a = make_agent()
    a.location.params.prep.type = ["Inj"]
    assert a.prep.last_dose_time is None

    model = make_model()

    rand_gen = FakeRandom(1.1)
    a.prep.enroll(rand_gen, 0)
    assert a.prep.last_dose_time == 0

    # make time pass
    model.time = 12
    a.prep.progress(model)
    assert a.prep.last_dose_time is None


@pytest.mark.unit
def test_cap_as_prob(make_agent, make_model, params):
    model = make_model(params)
    a = make_agent()
    a.location.params.prep.cap_as_prob = True

    model.run_random = FakeRandom(0.2)
    a.prep.initiate(model)
    assert not a.prep.active

    model.run_random = FakeRandom(0.0)
    a = make_agent()
    a.prep.initiate(model)
    assert a.prep.active

    model.run_random = FakeRandom(0.4)
    a.location.params.prep.target_model = ["Racial"]
    a.prep.active = False
    a.prep.initiate(model)
    assert not a.prep.active

    model.run_random = FakeRandom(0.3)
    a.prep.initiate(model)
    assert a.prep.active


@pytest.mark.unit
def test_get_acquisition_risk_multiplier(make_agent):
    a = make_agent()
    t = 0

    assert a.prep.get_acquisition_risk_multiplier(t, "sex") == 1.0

    a.prep.active = True
    a.prep.type = "Oral"
    a.prep.adherent = True
    a.prep.last_dose_time = -2

    assert a.prep.get_acquisition_risk_multiplier(t, "sex") == (
        1.0 - a.location.params.prep.efficacy.adherent
    )

    a.prep.type = "Inj"
    assert 0 < a.prep.get_acquisition_risk_multiplier(t, "sex") < 1.0
