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
    assert a.knowledge.active is False

    # remove all bonds compatible with agent 0. agent 0 fails
    for rel in copy(model.pop.relationships):
        if rel.bond_type in ["Inj", "SexInj"]:
            rel.unbond()
    with pytest.raises(ValueError) as excinfo:
        model.update_all_agents()
    assert "No agent zero!" in str(excinfo)


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
