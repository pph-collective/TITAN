import pytest

from titan.interactions import Injection
from titan.agent import Relationship

from conftest import FakeRandom


@pytest.mark.unit
def test_injection_transmission(make_model, make_agent):
    model = make_model()
    model.time = model.params.hiv.start_time + 2
    a = make_agent(race="white", DU="Inj", SO="HM")
    p = make_agent(race="white", DU="Inj", SO="HF")
    rel = Relationship(a, p, 10, bond_type="Inj")

    assert Injection.interact(model, rel) is False  # one agent must be HIV+

    a.hiv = True
    a.hiv_time = model.time - 1  # acute

    model.run_random = FakeRandom(-0.1)

    assert Injection.interact(model, rel)

    assert p.hiv


# TODO make this test the different bond types
