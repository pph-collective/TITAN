import pytest

from conftest import FakeRandom


@pytest.mark.unit
def test_external_exposure(make_model, make_agent, params):
    params.features.external_exposure = True

    model = make_model(params)
    external_exposure_agent = make_agent()
    external_exposure_agent.external_exposure.active = True
    model.pop.add_agent(external_exposure_agent)

    model.run_random = FakeRandom(-0.1)

    external_exposure_agent.external_exposure.update_agent(model)
    assert external_exposure_agent.hiv.active
