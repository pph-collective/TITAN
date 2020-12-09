import pytest

from copy import copy, deepcopy
import os
import math

from titan.model import *
from titan.agent import Relationship
from titan.features import HighRisk

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
    assert a.pca.awareness is False

    # remove all bonds compatible with agent 0. agent 0 fails
    for rel in copy(model.pop.relationships):
        if rel.bond_type in ["Inj", "SexInj"]:
            rel.unbond()
    with pytest.raises(ValueError) as excinfo:
        model.update_all_agents()
    assert "No agent zero!" in str(excinfo)


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

    a.incar.active = True
    assert model.agents_interact(rel) is False

    a.incar.active = False
    assert model.agents_interact(rel) is False  # neither HIV

    model.hiv_convert(a)
    model.hiv_convert(p)
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
    a.haart.active = True
    a.haart.adherent = True
    a.sex_role = "versatile"

    p = make_agent(race="white", SO="MSM")
    p.sex_role = "versatile"
    p.haart.active = True
    p.haart.adherent = True

    # test versatile-versatile relationship
    p_needle = (
        model.params.partnership.injection.transmission.haart_scaling.adherent
        * model.params.partnership.injection.transmission.base
    )
    p_sex = (
        model.params.partnership.sex.haart_scaling.MSM.adherent
        * model.params.partnership.sex.acquisition.MSM.versatile
    )
    scale = model.params.calibration.acquisition

    assert model.get_transmission_probability("injection", a, p) == p_needle * scale
    assert model.get_transmission_probability("sex", a, p) == p_sex * scale

    # test one vers agent, one receptive agent
    a.sex_role = "receptive"
    p_sex_ins = (
        model.params.partnership.sex.haart_scaling.MSM.adherent
        * model.params.partnership.sex.acquisition.MSM.insertive
    )
    p_sex_rec = (
        model.params.partnership.sex.haart_scaling.MSM.adherent
        * model.params.partnership.sex.acquisition.MSM.receptive
    )

    assert model.get_transmission_probability("sex", a, p) == p_sex_ins * scale
    assert model.get_transmission_probability("sex", p, a) == p_sex_rec * scale


@pytest.mark.unit
def test_hiv_convert(make_model, make_agent):
    model = make_model()
    a = make_agent()
    a.prep.active = True

    model.run_random = FakeRandom(-0.1)

    model.hiv_convert(a)

    assert a.hiv
    assert a.hiv_time == model.time
    assert a in model.pop.hiv_agents.members
    assert a.prep.active is False


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
    assert p.hiv_dx is False
    assert not p.partner_traced

    model.run_random = FakeRandom(-0.1)  # always less than param
    model.diagnose_hiv(a)

    assert a.hiv_dx
    assert a.hiv_dx_time == model.time
    assert p in a.get_partners()
    assert p.partner_traced
    assert p.trace_time == model.time

    assert p.hiv_dx is False
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
    a.haart.adherent = True  # .0051 prob

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
    a.haart.adherent = True  # .0051 prob

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
def test_die_and_replace_all(make_model, params):
    params.features.incar = False
    model = make_model(params)
    model.run_random = FakeRandom(0.0000001)  # always lower than death rate

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
    agent.incar.active = True
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
