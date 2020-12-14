import pytest
import os

from titan.location import *


@pytest.mark.unit
def test_location_init(params):
    location = "world"
    defn = params.classes.locations[location]

    world = Location(location, defn, params)

    assert world.ppl == 1.0
    assert world.name == location

    assert params.demographics.white.sex_type.WSW.drug_type.NonInj.hiv.aids.init == 0.1
    assert (
        world.params.demographics.white.sex_type.WSW.drug_type.NonInj.hiv.aids.init
        == 1.0
    )

    assert "white" in world.role_weights
    assert "black" in world.drug_weights
    assert "white" in world.pop_weights
    assert "values" in world.pop_weights["black"]
    assert "weights" in world.pop_weights["white"]

    assert len(world.edges) == 0


@pytest.mark.unit
def test_location_edge_init(params):
    defn = params.classes.locations["world"]
    location1 = Location("world", defn, params)
    location1.name = "world1"
    location2 = Location("world", defn, params)
    location2.name = "world2"

    edge = LocationEdge(location1, location2, 2.4)

    assert location1 in edge.edge
    assert location2 in edge.edge
    assert edge.distance == 2.4
    assert edge.id is not None

    location1a = Location("world", defn, params)
    location1a.name = "world1"
    with pytest.raises(AssertionError):
        LocationEdge(location1, location1a, 4.5)
