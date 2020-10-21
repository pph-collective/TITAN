import pytest

from titan.interactions import PCA
from titan.agent import Relationship

from conftest import FakeRandom

@pytest.mark.unit
def test_pca_interaction(make_model, make_agent):
    model = make_model()
    a = make_agent()
    p = make_agent()
    a.pca.opinion = 4
    p.pca.opinion = 2
    a.pca.awareness = True
    a.partners["SexInj"] = set()
    p.partners["SexInj"] = set()

    model.run_random = FakeRandom(1.0)

    model.pop.graph.add_edge(a, p)
    model.pop.graph.add_edge(a, "edge")

    model.time = 5

    rel = Relationship(a, p, 10, bond_type="SexInj")
    PCA.interact(model, rel)

    assert p.pca.awareness

    model.time += 1
    PCA.interact(model, rel)

    assert p.pca.opinion == 3
