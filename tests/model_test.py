import pytest

from copy import copy, deepcopy
import os
import math

from titan.model import *
from titan.agent import Agent, Relationship
from titan.parse_params import create_params

from conftest import FakeRandom


# ================================ MODEL TESTS =================================


@pytest.mark.unit
def test_model_init_error(params):
    params.model.seed.run = 0.5
    with pytest.raises(ValueError):
        HIVModel(params)


@pytest.mark.unit
def test_model_init(params):
    model = HIVModel(params)

    assert model.run_seed > 0
    assert model.pop.pop_seed > 0

    assert model.new_infections.num_members() == 0
    assert model.new_dx.num_members() == 0
    assert model.new_incar_release.num_members() == 0
    assert model.new_high_risk.num_members() == 0

    params.model.network.enable = False
    model = HIVModel(params)
    assert model.network_utils is None


@pytest.mark.unit
def test_update_AllAgents(make_model, make_agent):
    # make agent 0
    model = make_model()
    assert model.params.agent_zero.interaction_type == "injection"
    a = make_agent(race="white", DU="Inj")
    p = make_agent(race="white", DU="Inj")
    # make sure at least 1 relationship is compatible with agent 0 type
    Relationship(a, p, 10, bond_type="Inj")
    model.time = 1
    # update all agents passes with agent 0
    model.update_all_agents()
    assert a.prep_awareness is False

    # remove all bonds compatible with agent 0. agent 0 fails
    for rel in copy(model.pop.relationships):
        if rel.bond_type in ["Inj", "SexInj"]:
            rel.unbond()
    with pytest.raises(ValueError) as excinfo:
        model.update_all_agents()
    assert "No agent zero!" in str(excinfo)


@pytest.mark.unit
def test_pca_msmw(make_model, make_agent, params):

    # test update all agents for pca and msmw TODO separate tests
    model = make_model()
    a = make_agent(race="white", DU="Inj")
    p = make_agent(race="white", DU="Inj")
    model.pop.add_agent(a)
    model.pop.add_agent(p)
    msmw_agent = make_agent()
    msmw_agent.msmw = True
    model.pop.add_agent(msmw_agent)
    params.features.pca = True
    params.prep.pca.awareness.prob = 1.0

    model.time = 0
    model.update_all_agents()
    assert a.prep_awareness
    assert p.prep_awareness
    assert msmw_agent.hiv


@pytest.mark.unit
def test_initialize_random_trial_prep(make_model, params):
    params.features.prep = False
    params.vaccine.on_init = False
    model = make_model()
    model.run_random = FakeRandom(-0.1)
    model.time = 0
    model.initialize_random_trial()
    for agent in model.pop.all_agents:
        if not agent.hiv:
            assert agent.random_trial_enrolled
            assert agent.intervention_ever
            assert agent.prep


@pytest.mark.unit
def test_initialize_random_trial_pca_bridge(make_model, params):
    # pca trial
    model = make_model()
    params.features.pca = True
    model.initialize_random_trial()
    bridge_num = 0
    for comp in model.pop.connected_components():
        bridges = list(nx.bridges(comp))
        if bridges:
            bridge_num += 1
            for agent in comp.nodes():
                if agent.pca:
                    assert agent in [ag for ags in bridges for ag in ags]

    params.model.network.enable = False
    with pytest.raises(AssertionError) as excinfo:
        model.initialize_random_trial()
    assert "Network must be enabled for random trial" in str(excinfo)


@pytest.mark.unit
def test_agents_interact(make_model, make_agent):
    model = make_model()
    a = make_agent(race="white", SO="HM")
    p = make_agent(race="white", SO="HF")
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    model.run_random = FakeRandom(0.6)
    model.time = 0

    a.incar = True
    assert model.agents_interact(rel) is False

    a.incar = False
    assert model.agents_interact(rel) is False  # neither HIV

    a.hiv = True
    p.hiv = True
    assert model.agents_interact(rel) is False  # both HIV

    p.hiv = False

    assert model.agents_interact(rel)  # sex transmssion
    assert p.hiv is False  # but nothing happened (see skipped test)

    a.drug_type = "Inj"
    p.drug_type = "Inj"
    rel.bond_type = "Inj"

    model.run_random = FakeRandom(-0.1)

    assert model.agents_interact(rel)  # needle transmission
    assert p.hiv

    p.hiv = False
    model.run_random = FakeRandom(1.1)

    assert model.agents_interact(rel)  # needle and sex
    assert p.hiv is False  # but nothing happened


@pytest.mark.unit
def test_get_transmission_probability(make_model, make_agent):
    model = make_model()
    a = make_agent(race="white", SO="MSM")
    a.haart_adherence = 1
    a.sex_role = "versatile"

    p = make_agent(race="white", SO="MSM")
    p.sex_role = "versatile"
    p.haart_adherence = 1

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


@pytest.mark.unit
def test_needle_transmission(make_model, make_agent):
    model = make_model()
    a = make_agent(race="white", DU="Inj", SO="HM")
    p = make_agent(race="white", DU="Inj", SO="HF")

    with pytest.raises(AssertionError):
        model.injection_transmission(a, p)

    a.hiv = True
    a.hiv_time = 1  # acute

    model.run_random = FakeRandom(-0.1)

    model.injection_transmission(a, p)

    assert p.hiv


@pytest.mark.unit
def test_sex_transmission(make_model, make_agent):
    model = make_model()
    a = make_agent()
    a.sex_role = "insertive"
    p = make_agent()
    p.sex_role = "receptive"
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    a.hiv = True
    a.hiv_time = 1  # acute

    rel.total_sex_acts = 0
    model.params.calibration.acquisition = 10

    model.params.calibration.acquisition = 5
    model.params.calibration.sex.act = 10
    model.run_random = FakeRandom(0.6)

    # test partner becomes
    model.sex_transmission(rel)
    assert p.hiv


@pytest.mark.unit
def test_sex_transmission_do_nothing(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    with pytest.raises(ValueError):
        model.sex_transmission(rel)

    a.hiv = True
    p.hiv = True

    # test nothing happens
    model.sex_transmission(rel)


@pytest.mark.unit
def test_pca_interaction(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    a.prep_opinion = 4
    p.prep_opinion = 2
    a.prep_awareness = True
    a.partners["SexInj"] = set()
    p.partners["SexInj"] = set()

    model.run_random = FakeRandom(1.0)

    model.pop.graph.add_edge(a, p)
    model.pop.graph.add_edge(a, "edge")

    model.time = 5

    rel = Relationship(a, p, 10, bond_type="SexInj")
    model.pca_interaction(rel, force=True)

    assert p.prep_awareness

    model.time += 1
    model.pca_interaction(rel, force=True)

    assert p.prep_opinion == 3


@pytest.mark.unit
def test_hiv_convert(make_model, make_agent):
    model = make_model()
    a = make_agent()
    a.prep = True

    model.run_random = FakeRandom(-0.1)

    model.hiv_convert(a)

    assert a.hiv
    assert a.hiv_time == 1
    assert a in model.new_infections.members
    assert a in model.pop.hiv_agents.members
    assert a.prep is False


@pytest.mark.unit
def test_update_syringe_services(make_model):
    model = make_model()
    # make at least one agent PWID
    agent = next(iter(model.pop.all_agents))
    agent.drug_type = "Inj"
    if agent not in model.pop.pwid_agents.members:
        model.pop.pwid_agents.add_agent(agent)

    model.time = 3
    model.update_syringe_services()
    assert model.pop.pwid_agents

    for a in model.pop.all_agents:
        if a.drug_type == "Inj":
            assert a in model.pop.pwid_agents.members
            assert a.ssp


@pytest.mark.unit
def test_become_high_risk(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.become_high_risk(a, 10)

    assert a in model.pop.high_risk_agents.members
    assert a in model.new_high_risk.members
    assert a.high_risk
    assert a.high_risk_ever
    assert a.high_risk_time == 10

    model.params.features.high_risk = False
    assert not model.become_high_risk(a, 10)


@pytest.mark.unit
def test_update_high_risk(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # try to update when not high risk
    assert model.update_high_risk(a) is None

    a.high_risk = True
    a.high_risk_ever = True
    a.high_risk_time = 1
    model.pop.high_risk_agents.add_agent(a)

    model.update_high_risk(a)

    assert a.high_risk
    assert a.high_risk_time == 0
    assert a in model.pop.high_risk_agents

    model.update_high_risk(a)

    assert a.high_risk is False
    assert a.high_risk_ever is True
    assert a not in model.pop.high_risk_agents


@pytest.mark.unit
def test_incarcerate_not_hiv(make_model, make_agent):
    model = make_model()
    a = make_agent(SO="HM", race="white")
    p = make_agent(SO="HF", race="white")
    rel = Relationship(a, p, 10, "Sex")
    a.location.params.demographics.white.HM.incar.prob = 1.0

    model.incarcerate(a)
    assert a.incar


@pytest.mark.unit
def test_incarcerate_diagnosed(make_model, make_agent):
    model = make_model()
    model.time = 10
    a = make_agent(SO="HM", race="white")  # incarceration only for HM and HF?
    a.hiv = True
    a.hiv_dx = True
    a.partners["Sex"] = set()

    model.run_random = FakeRandom(-0.1)  # always less than params

    model.incarcerate(a)

    assert a.incar
    assert a.incar_time == 1
    assert a.haart
    assert a.haart_adherence == 5
    assert a.haart_time == 10


@pytest.mark.unit
def test_incarcerate_not_diagnosed(make_model, make_agent):
    model = make_model()
    a = make_agent(SO="HM", race="white")  # incarceration only for HM and HF?
    a.hiv = True
    a.partners["Sex"] = set()

    p = make_agent(SO="HF")
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    model.run_random = FakeRandom(-0.1)  # always less than params

    model.incarcerate(a)

    assert a.incar
    assert a.incar_time == 1
    assert a.hiv_dx

    assert p in model.pop.high_risk_agents.members
    assert p in model.new_high_risk.members
    assert p.high_risk
    assert p.high_risk_ever
    assert p.high_risk_time > 0


@pytest.mark.unit
def test_incarcerate_unincarcerate(make_model, make_agent):
    model = make_model()
    a = make_agent()
    a.target_partners = {bond: 0 for bond in model.params.classes.bond_types.keys()}
    a.mean_num_partners = copy(a.target_partners)

    a.incar = True
    a.incar_time = 2

    model.incarcerate(a)

    assert a.incar
    assert a.incar_time == 1

    model.incarcerate(a)

    assert a.incar is False
    assert a.incar_time == 0
    assert a in model.new_incar_release.members


@pytest.mark.unit
def test_diagnose_hiv(make_model, make_agent):
    model = make_model()
    model.params.partner_tracing.prob = 1.0
    model.time = 1
    a = make_agent()
    p = make_agent()
    p.hiv = True
    a.partners["Sex"].add(p)

    model.run_random = FakeRandom(1.1)  # always greater than param
    model.diagnose_hiv(a)

    assert a.hiv_dx is False
    assert a not in model.new_dx.members
    assert p.hiv_dx is False
    assert p not in model.new_dx.members
    assert not p.partner_traced

    model.run_random = FakeRandom(-0.1)  # always less than param
    model.diagnose_hiv(a)

    assert a.hiv_dx
    assert a in model.new_dx.members
    assert p in a.get_partners()
    assert p.partner_traced
    assert p.trace_time == model.time

    assert p.hiv_dx is False
    assert p not in model.new_dx.members
    model.params.demographics[p.race][p.sex_type].hiv.dx.prob = 0

    model.time = p.partner_traced + 1
    model.diagnose_hiv(p)
    assert p.hiv_dx
    assert p.partner_traced is False


@pytest.mark.unit
def test_diagnose_hiv_already_tested(make_model, make_agent):
    model = make_model()
    a = make_agent()

    a.hiv_dx = True

    model.run_random = FakeRandom(-0.1)  # always less than param
    model.diagnose_hiv(a)

    assert a.hiv_dx
    assert a not in model.new_dx.members


@pytest.mark.unit
def test_update_haart_t1(make_model, make_agent):
    model = make_model()
    model.time = 1
    a = make_agent(race="white")

    a.hiv = True

    # nothing happens, not tested
    model.update_haart(a)
    assert a.haart_adherence == 0
    assert a.haart is False

    # t0 agent initialized HAART
    a.hiv_dx = True

    # go on haart
    model.run_random = FakeRandom(
        -0.1
    )  # means this will always be less than params even though not physically possible in reality
    model.update_haart(a)

    assert a.haart_adherence == 5
    assert a.haart_time == 1
    assert a.haart

    # go off haart
    model.update_haart(a)

    assert a.haart_adherence == 0
    assert a.haart_time == 0
    assert a.haart is False


@pytest.mark.unit
def test_discontinue_prep_force(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # set up so the agent appears to be on PrEP
    a.prep = True
    num_prep = model.pop.prep_counts[a.race]
    model.pop.prep_counts[a.race] += 1

    model.discontinue_prep(a, True)

    assert a.prep is False
    assert num_prep == model.pop.prep_counts[a.race]


@pytest.mark.unit
def test_discontinue_prep_decrement_time(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.run_random = FakeRandom(1.1)

    # set up so the agent appears to be on PrEP
    a.prep = True
    num_prep = model.pop.prep_counts[a.race]
    model.pop.prep_counts[a.race] += 1

    model.discontinue_prep(a)

    assert a.prep
    assert num_prep + 1 == model.pop.prep_counts[a.race]


@pytest.mark.unit
def test_discontinue_prep_decrement_end(make_model, make_agent):
    model = make_model()
    a = make_agent(race="white")

    model.run_random = FakeRandom(-0.1)

    # set up so the agent appears to be on PrEP
    a.prep = True
    a.prep_type = "Oral"
    num_prep = model.pop.prep_counts[a.race]
    model.pop.prep_counts[a.race] += 1

    model.discontinue_prep(a)

    assert a.prep is False
    assert a.prep_type == ""
    assert num_prep == model.pop.prep_counts[a.race]


@pytest.mark.unit
def test_discontinue_prep_decrement_not_end(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.run_random = FakeRandom(1.1)

    # set up so the agent appears to be on PrEP
    a.prep = True
    a.prep_last_dose = 3
    a.prep_type = "Inj"
    num_prep = model.pop.prep_counts[a.race]
    model.pop.prep_counts[a.race] += 1

    model.discontinue_prep(a)

    assert a.prep
    assert a.prep_last_dose == 4  # 3 -> -1 -> +1 == 0 # Inj no longer in PrEP types
    assert num_prep + 1 == model.pop.prep_counts[a.race]
    assert a.prep_load > 0  # Inj no longer in PrEP types


@pytest.mark.unit
def test_discontinue_prep_decrement_inj_end(make_model, make_agent):
    model = make_model()
    a = make_agent()

    model.run_random = FakeRandom(1.1)

    # set up so the agent appears to be on PrEP
    a.prep = True
    a.prep_load = 0.4
    a.prep_last_dose = 12  # last step before hitting year mark and discontinuing
    a.prep_type = "Inj"
    num_prep = model.pop.prep_counts[a.race]
    model.pop.prep_counts[a.race] += 1

    model.discontinue_prep(a)

    assert a.prep is False
    assert a.prep_last_dose == 0  # 3 -> -1 -> +1 == 0 # Inj no longer in PrEP types
    assert num_prep == model.pop.prep_counts[a.race]
    assert a.prep_load == 0.0  # Inj no longer in PrEP types


@pytest.mark.unit
def test_initiate_prep_assertions(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # no PreP if already PreP
    a.prep = True
    assert model.initiate_prep(a) is None

    # no PrEP if already HIV
    a.prep = False
    a.hiv = True
    assert model.initiate_prep(a) is None


@pytest.mark.unit
def test_initiate_prep_force_adh(make_model, make_agent):
    model = make_model()
    a = make_agent()

    # forcing, adherant, inj
    model.run_random = FakeRandom(-0.1)
    model.initiate_prep(a, True)
    assert a.prep
    assert a in model.new_prep.members
    assert a.prep_adherence == 1
    # assert a.prep_load > 0.0 # Inj no longer in PrEP types
    assert a.prep_last_dose == 0


@pytest.mark.unit
def test_initiate_prep_force_non_adh(make_model, make_agent):
    model = make_model()
    a = make_agent()
    # forcing, non-adherant, inj
    model.run_random = FakeRandom(1.0)
    model.initiate_prep(a, True)
    assert a.prep
    assert a in model.new_prep.members
    assert a.prep_adherence == 0
    # assert a.prep_load > 0.0 # Inj no longer in PrEP types
    assert a.prep_last_dose == 0


@pytest.mark.unit
def test_initiate_prep_eligible(make_model, make_agent):
    model = make_model()

    # make sure there's room to add more prep agents
    a = make_agent(SO="HF")  # model is "CDCwomen"
    p = make_agent(DU="Inj")
    a.partners["Sex"] = set()
    p.partners["Sex"] = set()
    p.hiv_dx = True
    p.msmw = True
    model.time = 10
    a.location.params.prep.target = 1.0
    a.location.params.prep.target_model = ["cdc_women"]
    rel = Relationship(a, p, 10, bond_type="Sex")
    # non-forcing, adherant, inj
    model.run_random = FakeRandom(-0.1)
    model.initiate_prep(a)
    assert a.prep
    assert a in model.new_prep.members
    assert a.prep_adherence == 1
    # assert a.prep_load > 0.0 # Inj not in params prep_type anymore
    assert a.prep_last_dose == 0


@pytest.mark.unit
def test_progress_to_aids_error(make_agent, make_model):
    a = make_agent()
    a.hiv = False
    model = make_model()
    a.target_partners = {bond: 0 for bond in model.params.classes.bond_types.keys()}
    model.pop.add_agent(a)
    num_aids = sum([1 for agent in model.pop.hiv_agents if agent.aids])  # get baseline

    # test error case, agent must be HIV+
    with pytest.raises(AssertionError):
        model.progress_to_aids(a)

    assert sum([1 for agent in model.pop.hiv_agents if agent.aids]) == num_aids


@pytest.mark.unit
def test_progress_to_aids_nothing(make_agent, make_model):
    a = make_agent()
    a.hiv = True
    model = make_model()
    a.target_partners = {bond: 0 for bond in model.params.classes.bond_types.keys()}
    model.pop.add_agent(a)
    num_aids = sum([1 for agent in model.pop.hiv_agents if agent.aids])  # get baseline

    # test nothing case
    a.hiv = True
    a.haart_adherence = 1  # .0051 prob

    model.run_random = FakeRandom(0.9)  # no AIDS

    assert model.progress_to_aids(a) is None
    assert sum([1 for agent in model.pop.hiv_agents if agent.aids]) == num_aids
    assert a.aids is False


@pytest.mark.unit
def test_progress_to_aids_progress(make_agent, make_model):
    a = make_agent()
    a.hiv = True
    model = make_model()
    a.target_partners = {bond: 0 for bond in model.params.classes.bond_types.keys()}
    model.pop.add_agent(a)
    num_aids = sum([1 for agent in model.pop.hiv_agents if agent.aids])  # get baseline
    a.location.params.hiv.aids.prob = 1.0

    a.hiv = True
    a.haart_adherence = 1  # .0051 prob

    # test progress case
    model.run_random = FakeRandom(0.001)  # AIDS

    assert model.progress_to_aids(a) is None
    assert sum([1 for agent in model.pop.hiv_agents if agent.aids]) == num_aids + 1
    assert a.aids is True


@pytest.mark.unit
def test_die_and_replace_none(make_model):
    model = make_model()
    model.run_random = FakeRandom(0.999)  # always greater than death rate
    baseline_pop = copy(model.pop.all_agents.members)

    model.die_and_replace()

    ids = [a.id for a in baseline_pop]
    for agent in model.pop.all_agents.members:
        assert agent.id in ids


@pytest.mark.unit
def test_die_and_replace_all(make_model):
    model = make_model()
    model.run_random = FakeRandom(0.0000001)  # always lower than death rate

    # un-incarcerate everyone
    for agent in model.pop.all_agents.members:
        agent.incar = False

    baseline_pop = copy(model.pop.all_agents.members)
    old_ids = [a.id for a in baseline_pop]

    num_hm = len([x for x in baseline_pop if x.sex_type == "HM"])
    num_white = len([x for x in baseline_pop if x.race == "white"])
    num_pwid = len([x for x in baseline_pop if x.drug_type == "Inj"])

    model.die_and_replace()

    assert num_hm == len(
        [x for x in model.pop.all_agents.members if x.sex_type == "HM"]
    )
    assert num_white == len(
        [x for x in model.pop.all_agents.members if x.race == "white"]
    )
    assert num_pwid == len(
        [x for x in model.pop.all_agents.members if x.drug_type == "Inj"]
    )

    new_ids = [a.id for a in model.pop.all_agents.members]
    death_ids = [a.id for a in model.deaths]

    for agent in model.pop.all_agents.members:
        assert agent.id not in old_ids
        assert agent in model.pop.graph.nodes()

    for agent in baseline_pop:
        assert agent.id not in new_ids
        assert agent not in model.pop.graph.nodes()
        assert agent.id in death_ids


@pytest.mark.unit
def test_die_and_replace_incar(make_model):
    model = make_model()
    model.run_random = FakeRandom(0.0000001)  # always lower than death rate
    baseline_pop = copy(model.pop.all_agents.members)
    old_ids = [a.id for a in baseline_pop]

    agent = next(iter(model.pop.all_agents))
    agent.incar = True
    agent_id = agent.id

    model.die_and_replace()

    new_ids = [a.id for a in model.pop.all_agents.members]
    death_ids = [a.id for a in model.deaths]

    assert agent_id in old_ids
    assert agent_id not in death_ids
    assert agent_id in new_ids


@pytest.mark.unit
def test_timeline_scaling_default_def(make_model):
    model = make_model()
    original_params = deepcopy(model.params.__getstate__())
    model.time = 1
    model.timeline_scaling()

    assert original_params == model.params.__getstate__()


@pytest.mark.unit
def test_timeline_scaling_prep_def(make_model):
    model = make_model()
    scalar = 0.5
    model.params.timeline_scaling.timeline = ObjMap(
        {
            "scale": {
                "parameter": "prep|target",
                "start_time": 1,
                "stop_time": 3,
                "scalar": scalar,
            }
        }
    )
    original_prep_target = model.params.prep.target

    # scale the param
    model.time = 1
    model.timeline_scaling()

    assert math.isclose(
        original_prep_target * scalar, model.params.prep.target, abs_tol=0.001
    )

    # param still scaled
    model.time = 2
    model.timeline_scaling()

    assert math.isclose(
        original_prep_target * scalar, model.params.prep.target, abs_tol=0.001
    )

    # revert to original
    model.time = 3
    model.timeline_scaling()

    assert math.isclose(original_prep_target, model.params.prep.target, abs_tol=0.001)

    # still original
    model.time = 4
    model.timeline_scaling()

    assert math.isclose(original_prep_target, model.params.prep.target, abs_tol=0.001)
