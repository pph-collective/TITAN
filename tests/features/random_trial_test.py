import pytest
import networkx as nx  # type: ignore

from conftest import FakeRandom

from titan.features import RandomTrial


@pytest.mark.unit
def test_initialize_random_trial_prep_all(make_model, params):
    params.features.prep = True
    params.vaccine.on_init = False
    params.prep.target = 0
    params.random_trial.choice = "all"
    model = make_model(params)
    model.run_random = FakeRandom(-0.1)
    model.time = model.params.random_trial.start_time
    RandomTrial.update_pop(model)
    for agent in model.pop.all_agents:
        assert agent.random_trial.active
        if not agent.hiv.active:
            assert agent.random_trial.treated
            assert agent.prep.active


@pytest.mark.unit
def test_initialize_random_trial_prep_eigenvector(make_model, params):
    params.features.prep = True
    params.vaccine.on_init = False
    params.prep.target = 0
    params.hiv.start_time = 5
    params.random_trial.choice = "eigenvector"
    model = make_model(params)
    model.run_random = FakeRandom(-0.1)
    model.time = model.params.random_trial.start_time
    RandomTrial.update_pop(model)
    num_treated = 0
    num_suitable = 0
    for agent in model.pop.all_agents:
        assert agent.random_trial.active
        if agent.random_trial.treated:
            num_treated += 1
        if agent.random_trial.suitable:
            num_suitable += 1

    num_components = len(model.pop.connected_components())

    assert num_treated == num_components
    assert num_suitable == num_components


@pytest.mark.unit
def test_initialize_random_trial_prep_random(make_model, params):
    params.features.prep = True
    params.vaccine.on_init = False
    params.prep.target = 0
    params.hiv.start_time = 5
    params.random_trial.choice = "random"
    model = make_model(params)
    model.run_random = FakeRandom(-0.1)
    model.time = model.params.random_trial.start_time
    RandomTrial.update_pop(model)
    num_treated = 0
    num_suitable = 0
    for agent in model.pop.all_agents:
        assert agent.random_trial.active
        if agent.random_trial.treated:
            num_treated += 1
        if agent.random_trial.suitable:
            num_suitable += 1

    num_components = len(model.pop.connected_components())

    assert num_treated == num_components
    assert num_suitable == num_components


@pytest.mark.unit
def test_initialize_random_trial_pca_bridge(make_model, params):
    # knowledge bridge trial
    params.random_trial.treatment = "knowledge"
    model = make_model(params)
    model.time = model.params.random_trial.start_time
    RandomTrial.update_pop(model)
    bridge_num = 0
    components = model.pop.connected_components()
    assert len(components) > 0
    for comp in components:
        bridges = list(nx.bridges(comp))
        if bridges:
            bridge_num += 1
            for agent in comp.nodes():
                if agent.knowledge.active:
                    assert agent in [ag for ags in bridges for ag in ags]

    model.params.model.network.enable = False
    with pytest.raises(AssertionError) as excinfo:
        RandomTrial.update_pop(model)
    assert "Network must be enabled for random trial" in str(excinfo)
