import pytest
import networkx as nx  # type: ignore

from conftest import FakeRandom

from titan.features import RandomTrial


@pytest.mark.unit
def test_initialize_random_trial_prep(make_model, params):
    params.features.prep = False
    params.vaccine.on_init = False
    model = make_model()
    model.run_random = FakeRandom(-0.1)
    model.time = 0
    RandomTrial.update_pop(model)
    for agent in model.pop.all_agents:
        if not agent.hiv:
            assert agent.random_trial.active
            assert agent.random_trial.treated
            assert agent.prep.active


@pytest.mark.unit
def test_initialize_random_trial_pca_bridge(make_model, params):
    # pca trial
    model = make_model()
    params.features.pca = True
    RandomTrial.update_pop(model)
    bridge_num = 0
    for comp in model.pop.connected_components():
        bridges = list(nx.bridges(comp))
        if bridges:
            bridge_num += 1
            for agent in comp.nodes():
                if agent.pca.active:
                    assert agent in [ag for ags in bridges for ag in ags]

    params.model.network.enable = False
    with pytest.raises(AssertionError) as excinfo:
        RandomTrial.update_pop(model)
    assert "Network must be enabled for random trial" in str(excinfo)
