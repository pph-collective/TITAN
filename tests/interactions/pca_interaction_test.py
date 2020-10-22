import pytest

from titan.interactions.pca import *
from titan.agent import Relationship
from titan import utils

from conftest import FakeRandom


@pytest.mark.unit
def test_pca_interaction(make_model, make_agent):
    model = make_model()
    model.params.features.pca = True
    a = make_agent()
    p = make_agent()
    a.pca.opinion = 4
    p.pca.opinion = 2
    a.pca.awareness = True
    a.partners["SexInj"] = set()
    p.partners["SexInj"] = set()
    rel = Relationship(a, p, 10, bond_type="SexInj")

    model.run_random = FakeRandom(-0.1)

    model.pop.graph.add_edge(a, p)
    model.pop.graph.add_edge(a, "edge")

    model.time = 5

    # make partner aware via dissemination
    PCA.interact(model, rel)

    assert p.pca.awareness

    # influence partner's opinion
    model.time += 1
    PCA.interact(model, rel)

    assert p.pca.opinion == 3


@pytest.mark.unit
def test_pca_knowledge_transmission_probability(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    rel = Relationship(a, p, 10, bond_type="Social")

    # knowledge transmission
    assert knowledge_transmission_probability(model, rel, 0) == 0
    assert (
        knowledge_transmission_probability(model, rel, 1)
        == model.params.pca.knowledge.transmission
    )
    assert knowledge_transmission_probability(model, rel, 2) == 1.0 - utils.binom_0(
        2, model.params.pca.knowledge.transmission
    )

    # opinion transmission
    a.pca.awareness = True
    p.pca.awareness = True

    assert knowledge_transmission_probability(model, rel, 0) == 0
    assert (
        knowledge_transmission_probability(model, rel, 1)
        == model.params.pca.opinion.transmission
    )
    assert knowledge_transmission_probability(model, rel, 2) == 1.0 - utils.binom_0(
        2, model.params.pca.opinion.transmission
    )


@pytest.mark.unit
def test_pca_knowledge_dissemination(make_model, make_agent):
    model = make_model()
    a = make_agent()

    a.pca.opinion = 0
    model.run_random = FakeRandom(1)
    knowledge_dissemination(model, a)

    assert a.pca.awareness
    assert not a.prep.active

    a.pca.awareness = False
    a.pca.opinion = 5
    model.run_random = FakeRandom(-0.1)
    knowledge_dissemination(model, a)

    assert a.pca.awareness
    assert a.prep.active


@pytest.mark.unit
def test_pca_influence(make_model, make_agent):
    model = make_model()
    model.run_random = FakeRandom(-0.1)
    rel = next(iter(model.pop.relationships))  # get a relationship from the pop

    rel.agent1.prep.active = False
    rel.agent2.prep.active = False
    rel.agent1.hiv = False
    rel.agent2.hiv = False

    rel.agent1.pca.opinion = 2
    rel.agent2.pca.opinion = 5

    influence(model, rel)

    # did agent1's opinion change? otherwise, try the other way
    if rel.agent1.pca.opinion > 2:
        assert rel.agent1.prep.active
    else:
        rel.agent1.pca.opinion = 5
        rel.agent2.pca.opinion = 2

        influence(model, rel)

        assert rel.agent2.prep.active
