import pytest

from conftest import FakeRandom

from titan.agent import Relationship
from titan import utils
from titan.exposures import influence


@pytest.mark.unit
def test_knolwedge_transmission_probability(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    rel = Relationship(a, p, 10, bond_type="Social")

    # knowledge transmission
    assert a.knowledge.get_transmission_probability(model, "pca", p, 0) == 0
    assert (
        a.knowledge.get_transmission_probability(model, "pca", p, 1)
        == model.params.knowledge.prob
    )
    assert a.knowledge.get_transmission_probability(
        model, "pca", p, 2
    ) == 1.0 - utils.binom_0(2, model.params.knowledge.prob)

    # opinion transmission
    a.knowledge.active = True
    p.knowledge.active = True

    assert a.knowledge.get_transmission_probability(model, "pca", p, 0) == 0
    assert (
        a.knowledge.get_transmission_probability(model, "pca", p, 1)
        == model.params.knowledge.opinion.prob
    )
    assert a.knowledge.get_transmission_probability(
        model, "pca", p, 2
    ) == 1.0 - utils.binom_0(2, model.params.knowledge.opinion.prob)


@pytest.mark.unit
def test_knowledge_convert(make_model, make_agent):
    model = make_model()
    a = make_agent()

    a.knowledge.opinion = 0
    model.run_random = FakeRandom(1)
    a.knowledge.convert(model)

    assert a.knowledge.active
    assert not a.prep.active

    a.knowledge.active = False
    a.knowledge.opinion = 5
    model.run_random = FakeRandom(-0.1)
    a.knowledge.convert(model)

    assert a.knowledge.active
    assert a.prep.active


@pytest.mark.unit
def test_knowledge_influence(make_model, make_agent):
    model = make_model()
    model.run_random = FakeRandom(-0.1)
    rel = next(iter(model.pop.relationships))  # get a relationship from the pop

    rel.agent1.prep.active = False
    rel.agent2.prep.active = False
    rel.agent1.hiv.active = False
    rel.agent2.hiv.active = False

    rel.agent1.knowledge.opinion = 2
    rel.agent2.knowledge.opinion = 5

    influence(model, rel)

    # did agent1's opinion change? otherwise, try the other way
    if rel.agent1.knowledge.opinion > 2:
        assert rel.agent1.prep.active
    else:
        rel.agent1.knowledge.opinion = 5
        rel.agent2.knowledge.opinion = 2

        influence(model, rel)

        assert rel.agent2.prep.active


@pytest.mark.unit
def test_knowledge_exposure(make_model, make_agent, params):

    # test update all agents for pca and msmw
    params.exposures.knowledge = True
    params.knowledge.prob = 0.0

    model = make_model(params)
    a = make_agent(race="white", DU="Inj")
    a.knowledge.opinion = 4
    model.pop.add_agent(a)

    assert not a.knowledge.active
    assert not a.prep.active

    model.time = model.params.knowledge.start_time
    model.run_random = FakeRandom(-0.1)
    a.knowledge.update_agent(model)

    assert a.knowledge.active
    assert a.prep.active
