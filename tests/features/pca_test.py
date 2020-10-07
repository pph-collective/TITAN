import pytest

from conftest import FakeRandom


@pytest.mark.unit
def test_pca(make_model, make_agent, params):

    # test update all agents for pca and msmw TODO separate tests
    params.features.pca = True
    params.pca.awareness.prob = 0.0

    model = make_model(params)
    a = make_agent(race="white", DU="Inj")
    model.pop.add_agent(a)

    assert not a.pca.awareness
    assert not a.prep.active

    model.time = model.params.pca.start_time
    model.run_random = FakeRandom(-0.1)
    a.pca.update_agent(model)

    assert a.pca.awareness
    assert a.prep.active
