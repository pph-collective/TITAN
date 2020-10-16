import pytest

from titan.features import SyringeServices


@pytest.mark.unit
def test_update_syringe_services(make_model):
    model = make_model()
    # make at least one agent PWID
    agent = next(iter(model.pop.all_agents))
    agent.drug_type = "Inj"
    if agent not in model.pop.pwid_agents.members:
        model.pop.pwid_agents.add_agent(agent)

    model.time = 3
    SyringeServices.update_pop(model)
    assert model.pop.pwid_agents

    for a in model.pop.all_agents:
        if a.drug_type == "Inj":
            assert a in model.pop.pwid_agents.members
            assert a.syringe_services.active
