import pytest

from conftest import FakeRandom


@pytest.mark.unit
def test_msmw(make_model, make_agent, params):
    params.features.msmw = True

    model = make_model(params)
    msmw_agent = make_agent()
    msmw_agent.msmw.active = True
    model.pop.add_agent(msmw_agent)

    model.run_random = FakeRandom(-0.1)

    msmw_agent.msmw.update_agent(model)
    assert msmw_agent.hiv
