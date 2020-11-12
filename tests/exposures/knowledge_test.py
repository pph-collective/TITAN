import pytest


@pytest.mark.skip
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
    a.knowledge.active = True
    p.knowledge.active = True

    assert knowledge_transmission_probability(model, rel, 0) == 0
    assert (
        knowledge_transmission_probability(model, rel, 1)
        == model.params.knowledge.opinion.transmission
    )
    assert knowledge_transmission_probability(model, rel, 2) == 1.0 - utils.binom_0(
        2, model.params.knowledge.opinion.transmission
    )


@pytest.mark.skip
def test_pca_knowledge_dissemination(make_model, make_agent):
    model = make_model()
    a = make_agent()

    a.knowledge.opinion = 0
    model.run_random = FakeRandom(1)
    knowledge_dissemination(model, a)

    assert a.knowledge.active
    assert not a.prep.active

    a.knowledge.active = False
    a.knowledge.opinion = 5
    model.run_random = FakeRandom(-0.1)
    knowledge_dissemination(model, a)

    assert a.knowledge.active
    assert a.prep.active


@pytest.mark.skip
def test_pca_influence(make_model, make_agent):
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


@pytest.mark.skip
def test_pca(make_model, make_agent, params):

    # test update all agents for pca and msmw TODO separate tests
    params.features.pca = True
    params.knowledge.active.prob = 0.0

    model = make_model(params)
    a = make_agent(race="white", DU="Inj")
    model.pop.add_agent(a)

    assert not a.knowledge.active
    assert not a.prep.active

    model.time = model.params.pca.start_time
    model.run_random = FakeRandom(-0.1)
    a.pca.update_agent(model)

    assert a.knowledge.active
    assert a.prep.active
