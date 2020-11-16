import pytest
from copy import copy

from conftest import FakeRandom

from titan.features import Incar, HighRisk
from titan.agent import Relationship


@pytest.mark.unit
def test_incarcerate_unincarcerate(make_model, make_agent):
    model = make_model()
    a = make_agent()
    a.target_partners = {bond: 0 for bond in model.params.classes.bond_types.keys()}
    a.mean_num_partners = copy(a.target_partners)

    a.incar.active = True
    a.incar.duration = 2

    a.incar.update_agent(model)

    assert a.incar.active
    assert a.incar.duration == 1

    a.incar.update_agent(model)

    assert a.incar.active is False
    assert a.incar.duration == 0
    assert a in Incar.new_releases


@pytest.mark.unit
def test_incarcerate_not_diagnosed(make_model, make_agent):
    model = make_model()
    a = make_agent(SO="HM", race="white")  # incarceration only for HM and HF?
    a.hiv = True
    a.partners["Sex"] = set()

    p = make_agent(SO="HF")
    p.partners["Sex"] = set()
    rel = Relationship(a, p, 10, bond_type="Sex")

    model.run_random = FakeRandom(-0.1)  # always less than params

    a.incar.update_agent(model)

    assert a.incar.active
    assert a.incar.duration == 1
    assert a.hiv_dx

    assert p.high_risk.active
    assert p.high_risk.ever
    assert p.high_risk.duration > 0
    assert p.high_risk.time == model.time


@pytest.mark.unit
def test_incarcerate_not_hiv(make_model, make_agent):
    model = make_model()
    a = make_agent(SO="HM", race="white")
    p = make_agent(SO="HF", race="white")
    rel = Relationship(a, p, 10, "Sex")
    a.location.params.demographics.white.HM.incar.prob = 1.0

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
    assert a.incar.duration == 1
    assert a.haart.active
    assert a.haart.adherent is True
