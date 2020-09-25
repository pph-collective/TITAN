import pytest


@pytest.mark.unit
def test_pca_msmw(make_model, make_agent, params):

    # test update all agents for pca and msmw TODO separate tests
    model = make_model()
    a = make_agent(race="white", DU="Inj")
    p = make_agent(race="white", DU="Inj")
    model.pop.add_agent(a)
    model.pop.add_agent(p)
    msmw_agent = make_agent()
    msmw_agent.msmw.active = True
    model.pop.add_agent(msmw_agent)
    params.features.pca = True
    params.pca.awareness.prob = 1.0

    model.time = 0
    model.update_all_agents()
    assert a.pca.awareness
    assert p.pca.awareness
    assert msmw_agent.hiv
