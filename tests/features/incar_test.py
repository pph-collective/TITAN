import pytest
from copy import copy

from conftest import FakeRandom

from titan.features import Incar, HighRisk
from titan.agent import Relationship


@pytest.mark.unit
def test_incarcerate_unincarcerate(make_model, make_agent):
    model = make_model()
    model.run_random = FakeRandom(-0.1)
    a = make_agent()
    a.target_partners = {bond: 0 for bond in model.params.classes.bond_types.keys()}
    a.mean_num_partners = copy(a.target_partners)

    a.incar.active = True
    a.incar.release_time = model.time + 2
    a.hiv = True
    a.haart.active = True

    model.time += 1
    a.incar.update_agent(model)

    assert a.incar.active

    model.time += 1
    a.incar.update_agent(model)

    assert a.incar.active is False
    assert a.haart.active is False
    assert a.haart.adherent is False


@pytest.mark.unit
def test_incarcerate_not_diagnosed(make_model, make_agent):
    model = make_model()
    a = make_agent(SO="HM", race="white")  # incarceration only for HM and HF?
    a.hiv = True

    model.run_random = FakeRandom(-0.1)  # always less than params

    a.incar.update_agent(model)

    assert a.incar.active
    assert a.incar.release_time == model.time + 1
    assert a.hiv_dx


@pytest.mark.unit
def test_incarcerate_not_hiv(make_model, make_agent):
    model = make_model()
    a = make_agent(SO="HM", race="white")
    p = make_agent(SO="HF", race="white")
    rel = Relationship(a, p, 10, "Sex")
    a.location.params.demographics.white.sex_type.HM.incar.prob = 1.0

    a.incar.update_agent(model)
    assert a.incar.active


@pytest.mark.unit
def test_incarcerate_diagnosed(make_model, make_agent):
    model = make_model()
    model.time = 10
    a = make_agent(SO="HM", race="white")  # incarceration only for HM and HF?
    a.hiv = True
    a.hiv_dx = True
    a.partners["Sex"] = set()

    model.run_random = FakeRandom(-0.1)  # always less than params

    a.incar.update_agent(model)

    assert a.incar.active
    assert a.incar.release_time == model.time + 1
    assert a.haart.active
    assert a.haart.adherent is True

    # Goes on haart but nonadherent
    a = make_agent(SO="HM", race="white")
    a.location.params.incar.haart.adherence = -1.0
    a.hiv = True
    a.hiv_dx = True
    a.partners["Sex"] = set()
    model.run_random = FakeRandom(-0.1)  # between haart adherence and other params
    a.incar.update_agent(model)
    assert a.incar.active
    assert a.incar.release_time == model.time + 1
    assert a.haart.active
    assert not a.haart.adherent
