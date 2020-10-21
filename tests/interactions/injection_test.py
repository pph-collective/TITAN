import pytest

from titan.interactions import Injection
from titan.agent import Relationship

from conftest import FakeRandom

@pytest.mark.unit
def test_injection_transmission(make_model, make_agent):
    model = make_model()
    a = make_agent(race="white", DU="Inj", SO="HM")
    p = make_agent(race="white", DU="Inj", SO="HF")
    rel = Relationship(a, p, 10, bond_type="Inj")

    with pytest.raises(AssertionError):
        Injection.interact(model, rel)

    a.hiv = True
    a.hiv_time = model.time - 1  # acute

    model.run_random = FakeRandom(-0.1)

    Injection.interact(model, rel)

    assert p.hiv
