import pytest

from titan.interactions.pca import *
from titan.agent import Relationship
from titan import utils

from conftest import FakeRandom


@pytest.mark.unit
def test_pca_interaction(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    a.knowledge.opinion = 4
    p.knowledge.opinion = 2
    a.knowledge.active = True
    a.partners["SexInj"] = set()
    p.partners["SexInj"] = set()
    rel = Relationship(a, p, 10, bond_type="SexInj")

    model.run_random = FakeRandom(-0.1)

    model.pop.graph.add_edge(a, p)
    model.pop.graph.add_edge(a, "edge")

    model.time = 5

    # make partner aware via dissemination
    PCA.interact(model, rel)

    assert p.knowledge.active

    # influence partner's opinion
    model.time += 1
    PCA.interact(model, rel)

    assert p.knowledge.opinion == 3


@pytest.mark.unit
def test_pca_num_acts(make_model, make_agent):
    model = make_model()
    model.run_random = FakeRandom(-0.1)
    a = make_agent()
    p = make_agent()
    rel = Relationship(a, p, 10, bond_type="Social")

    assert (
        PCA.get_num_acts(model, rel)
        == model.params.partnership.pca.frequency.Social[1].min
    )
