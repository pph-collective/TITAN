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
        TITAN(params)


@pytest.mark.unit
def test_model_init(params):
    model = TITAN(params)

    assert model.run_seed > 0
    assert model.pop.pop_seed > 0


@pytest.mark.unit
def test_update_all_agents(make_model, make_agent):
    # make agent 0
    model = make_model()
    assert model.params.agent_zero.interaction_type == "injection"
    model.time = 1
    model.params.features.exit_enter = False
    # update all agents passes with agent 0
    model.update_all_agents()

    # remove all bonds compatible with agent 0. agent 0 fails
    for rel in copy(model.pop.relationships):
        if rel.bond_type in ["Inj", "SexInj"]:
            rel.unbond()
            model.pop.remove_relationship(rel)

    # no target partners available
    for agent in model.pop.all_agents:
        for k in agent.target_partners.keys():
            agent.target_partners[k] = 0

    model.reset_trackers()

    with pytest.raises(ValueError) as excinfo:
        model.update_all_agents()

    assert "No agent zero!" in str(excinfo)


@pytest.mark.unit
def test_death_none(make_model):
    model = make_model()
    model.run_random = FakeRandom(0.999)
    baseline_pop = copy(model.pop.all_agents.members)
    model.exit()

    ids = [a.id for a in baseline_pop]
    for agent in model.pop.all_agents.members:
        assert agent.id in ids


@pytest.mark.unit
def test_death_incar(make_model, params):
    params.demographics.white.sex_type.MSM.incar.init = 1.0
    model = make_model(params)
    model.run_random = FakeRandom(0.00000001)
    incar_pop = {agent for agent in model.pop.all_agents if agent.incar.active}
    non_incar = {agent for agent in model.pop.all_agents if not agent.incar.active}
    model.exit()

    # no incarcerated agents exit if incar is ignored
    assert incar_pop == model.pop.all_agents.members
    # check that list of removed agents is updated
    for agent in non_incar:
        assert agent in model.exits["death"]

    model.reset_trackers()
    model.params.classes.exit.death.ignore_incar = False
    model.exit()
    assert len(model.pop.all_agents.members) == 0


@pytest.mark.unit
def test_death_all(make_model, params):
    params.features.incar = False
    model = make_model()
    model.run_random = FakeRandom(0.0000000001)
    model.exit()

    assert len(model.pop.all_agents.members) == 0


@pytest.mark.unit
def test_dropout_none(make_model, params):
    params.exit_enter = ObjMap(
        {"dropout": {"exit_class": "migrate", "entry_class": "none"}}
    )
    model = make_model(params)
    init_ag = {agent for agent in model.pop.all_agents.members}
    # no dropouts
    model.run_random = FakeRandom(1.0)
    model.exit()

    end_ag = {agent for agent in model.pop.all_agents.members}
    assert init_ag == end_ag


@pytest.mark.unit
def test_dropout_all(make_model, params):
    params.exit_enter = ObjMap(
        {"dropout": {"exit_class": "migrate", "entry_class": "none"}}
    )
    model = make_model(params)
    init_ppl = copy(model.pop.all_agents.members)
    model.run_random = FakeRandom(0.000000001)
    model.exit()

    for agent in init_ppl:
        assert agent in model.exits["migrate"]
    for exit_type, agents in model.exits.items():
        if exit_type != "migrate":
            assert agents == []


@pytest.mark.unit
def test_ageout(make_model, params):
    params.exit_enter = ObjMap(
        {"dropout": {"exit_class": "age_out", "entry_class": "none"}}
    )
    model = make_model(params)
    init_ppl = copy(model.pop.all_agents.members)
    model.exit()

    assert {agent for agent in model.pop.all_agents.members} == init_ppl

    model.params.classes.exit.age_out.age = 45
    model.exit()

    for agent in model.pop.all_agents.members:
        assert agent.age <= 45
    for agent in model.exits["age_out"]:
        assert agent.age > 45

    model.reset_trackers()
    model.params.classes.exit.age_out.age = 0
    model.exit()
    assert not len(model.pop.all_agents.members)


@pytest.mark.unit
def test_exit_none(make_model, params):
    params.exit_enter.death.exit_class = "none"
    params.features.incar = False
    model = make_model(params)
    init_ppl = copy(model.pop.all_agents.members)
    model.run_random = FakeRandom(0.000000001)

    assert not model.exit()
    assert model.pop.all_agents.members == init_ppl


@pytest.mark.unit
def test_new_agent_no_exit(make_model, params):
    params.exit_enter.death.exit_class = "none"
    params.exit_enter.death.entry_class = "new_ag"
    model = make_model(params)
    model.run_random = FakeRandom(0.5)
    init_ppl = copy(model.pop.all_agents.members)

    model.enter()
    assert len(init_ppl) * 2 == len(model.pop.all_agents.members)

    init_black = 0
    init_white = 0
    new_black = 0
    new_white = 0
    for agent in init_ppl:
        if agent.race == "black":
            init_black += 1
        elif agent.race == "white":
            init_white += 1
    for agent in model.pop.all_agents.members:
        if agent.race == "black":
            new_black += 1
        if agent.race == "white":
            new_white += 1
    # check that ratios are the same
    assert init_black / len(init_ppl) == new_black / len(model.pop.all_agents.members)
    assert init_white / len(init_ppl) == new_white / len(model.pop.all_agents.members)


@pytest.mark.unit
def test_new_agent(make_model, params):
    params.exit_enter.death.exit_class = "death"
    params.exit_enter.death.entry_class = "new_ag"
    model = make_model(params)
    model.run_random = FakeRandom(0.0000000001)
    init_ppl = copy(model.pop.all_agents.members)

    model.exit()
    assert len(init_ppl) > model.pop.all_agents.num_members()
    model.enter()
    assert len(init_ppl) == model.pop.all_agents.num_members()

    model.reset_trackers()
    model.params.classes.enter.new_ag.age_in = True
    model.exit()
    model.enter()
    for agent in model.pop.all_agents.members:
        if not agent.incar.active:
            assert agent.age == params.classes.enter.new_ag.age


@pytest.mark.unit
def test_replace(make_model, params):
    params.features.incar = False
    model = make_model()
    init_ppl = copy(model.pop.all_agents.members)
    # run with no replace list
    model.run_random = FakeRandom(0.0000001)
    model.enter()
    assert model.pop.all_agents.members == init_ppl

    # all agents exit and replaced
    model.exit()
    model.enter()
    assert model.pop.all_agents.members != init_ppl
    assert model.pop.all_agents.num_members() == len(init_ppl)

    model.reset_trackers()
    # all agents exit but none replaced
    model.exit()
    model.run_random = FakeRandom(1.0)
    model.enter()
    assert not model.pop.all_agents.members


@pytest.mark.unit
def test_agein(make_model, params):
    params.exit_enter.death.entry_class = "age_in"
    params.features.incar = False
    model = make_model(params)
    init_ppl = copy(model.pop.all_agents.members)
    model.run_random = FakeRandom(0.00000001)
    model.exit()

    model.enter()
    assert model.pop.all_agents.members != init_ppl
    assert model.pop.all_agents.num_members() == len(init_ppl)
    for agent in model.pop.all_agents.members:
        assert agent.age == 16


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
                "parameter": "prep|cap",
                "start_time": 1,
                "stop_time": 3,
                "scalar": scalar,
            }
        }
    )
    original_prep_target = model.params.prep.cap

    # scale the param
    model.time = 1
    model.timeline_scaling()

    assert math.isclose(
        original_prep_target * scalar, model.params.prep.cap, abs_tol=0.001
    )

    # param still scaled
    model.time = 2
    model.timeline_scaling()

    assert math.isclose(
        original_prep_target * scalar, model.params.prep.cap, abs_tol=0.001
    )

    # revert to original
    model.time = 3
    model.timeline_scaling()

    assert math.isclose(original_prep_target, model.params.prep.cap, abs_tol=0.001)

    # still original
    model.time = 4
    model.timeline_scaling()

    assert math.isclose(original_prep_target, model.params.prep.cap, abs_tol=0.001)
