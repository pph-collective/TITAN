import pytest

from titan.interactions import Injection
from titan.agent import Relationship
from titan.features import SyringeServices

from conftest import FakeRandom


@pytest.mark.unit
def test_injection_transmission(make_model, make_agent):
    model = make_model()
    model.np_random = FakeRandom(1.0)
    model.run_random = FakeRandom(-0.1)
    model.time = model.params.hiv.start_time + 2
    a = make_agent(race="black", DU="Inj", SO="HM")
    p = make_agent(race="black", DU="Inj", SO="HF")
    rel = Relationship(a, p, 10, bond_type="Inj")

    Injection.interact(model, rel)
    assert a.hiv.active is False
    assert p.hiv.active is False

    a.hiv.active = True
    a.hiv.time = model.time - 1  # acute

    Injection.interact(model, rel)

    assert p.hiv.active


@pytest.mark.unit
def test_injection_num_acts(make_model, make_agent):
    model = make_model()
    model.np_random = FakeRandom(0.5)
    a = make_agent()
    p = make_agent()
    a.drug_type = "Inj"
    p.drug_type = "Inj"
    rel = Relationship(a, p, 10, bond_type="Inj")

    # set to a high number to ensure above zero
    a.location.params.demographics[a.race].sex_type[a.sex_type].injection.num_acts = 100
    p.location.params.demographics[p.race].sex_type[p.sex_type].injection.num_acts = 200
    assert Injection.get_num_acts(model, rel) > 0

    a.syringe_services.active = True
    SyringeServices.enrolled_risk = 0.0

    assert Injection.get_num_acts(model, rel) == 0

    a.syringe_services.active = False
    a.hiv.active = True
    a.hiv.dx = True
    model.params.hiv.dx.risk_reduction.injection = 1.0

    assert Injection.get_num_acts(model, rel) == 0
    assert p.hiv


@pytest.mark.unit
def test_injection_get_num_acts_do_nothing(make_model, make_agent):
    model = make_model()
    model.time = model.params.hiv.start_time + 2
    a = make_agent(race="white", DU="Inj", SO="HM")
    p_inj = make_agent(race="white", DU="Inj", SO="HF")
    rel_Inj = Relationship(a, p_inj, 10, bond_type="Inj")

    model.run_random = FakeRandom(-0.1)

    assert Injection.get_num_acts(model, rel_Inj) > 0

    model.run_random = FakeRandom(1.1)
    assert Injection.get_num_acts(model, rel_Inj) == 0
