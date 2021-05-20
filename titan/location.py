from typing import Optional, Set, Dict, List, Any
from copy import deepcopy
import math
import os
import csv

from .parse_params import ObjMap
from . import utils


class Location:
    def __init__(self, name: str, defn: ObjMap, params: ObjMap):
        """
        This class constructs and represents a location within the model.  A location
            can have an arbitrary geographic granularity.

        args:
            name: name of the location
            defn: definition for this location
            params: model parameters
        """
        # location properties
        self.name = name
        self.params = self.create_params(params)
        self.ppl = defn.ppl  # percent of overall population assigned to this location
        self.category = defn.category  # arbitrary category, can be used for migration

        # value/weight maps needed for creating new agents in this location
        self.pop_weights: Dict[str, Dict[str, List[Any]]] = {}
        self.role_weights: Dict[str, Dict] = {}
        self.drug_weights: Dict[str, Dict] = {}
        self.init_weights()

        self.migration_weights: Dict[str, Any] = {}

        self.neighbors: Set[str] = set()  # or maybe edges instead

    def __str__(self):
        return self.name

    def __repr__(self):
        return f"'{self.name}'"

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return self.name != other.name

    def __hash__(self):
        return hash(self.name)

    def create_params(self, params: ObjMap) -> ObjMap:
        """
        Scale or override the generic parameters with any location based scaling from params.location.scaling

        args:
            params: model parameters

        returns:
            new parameter object with scaled values for this location
        """
        new_params = deepcopy(params)

        defns = new_params.location.scaling[self.name]
        for param_path, defn in defns.items():
            if param_path != "ls_default":
                if defn.field == "scalar":
                    utils.scale_param(new_params, param_path, defn.scalar)
                elif defn.field == "override":
                    utils.override_param(new_params, param_path, defn.override)

        return new_params

    def init_weights(self):
        """
        Create the containers to hold values and weights for randomly selecting:

        * sex_role
        * drug_type
        * race
        * sex_type
        """

        def init_weight_dict(d, item):
            d[item] = {}
            d[item]["values"] = []
            d[item]["weights"] = []

        def add_weight(d, v, w):
            d["values"].append(v)
            d["weights"].append(w)

        total_ppl = 0
        for race, race_param in self.params.demographics.items():
            self.role_weights[race] = {}
            self.drug_weights[race] = {}
            init_weight_dict(self.pop_weights, race)
            total_ppl += race_param.ppl
            for st, st_param in race_param.sex_type.items():
                add_weight(self.pop_weights[race], st, st_param.ppl)
                init_weight_dict(self.role_weights[race], st)
                init_weight_dict(self.drug_weights[race], st)
                for role, prob in st_param.sex_role.init.items():
                    add_weight(self.role_weights[race][st], role, prob)
                for dt, dt_param in st_param.drug_type.items():
                    add_weight(self.drug_weights[race][st], dt, dt_param.ppl)

                assert math.isclose(
                    sum(self.role_weights[race][st]["weights"]), 1, abs_tol=0.001
                ), f"{self.name}'s' {race} {st} role weights must add to 1"
                assert math.isclose(
                    sum(self.drug_weights[race][st]["weights"]), 1, abs_tol=0.001
                ), f"ppl of {self.name}'s' {race} {st} drug_types must add to 1"

            assert math.isclose(
                sum(self.pop_weights[race]["weights"]), 1, abs_tol=0.001
            ), f"ppl of {self.name}'s' {race} sex_types must add to 1"

        assert math.isclose(
            total_ppl, 1, abs_tol=0.001
        ), f"ppl of {self.name}'s' races must add to 1"


# LocationEdges are very much a WIP and not actually used anywhere yet
# outstanding questions:
# * should edges be directed? bi-drectional? uni-drectional, but both directions housed in the same edge?
# * what attributes do edges need?
# * how will mobility be implemented?
# * assorting?
class LocationEdge:

    next_edge_id = 0

    @classmethod
    def update_id_counter(cls, last_id: int):
        cls.next_edge_id = last_id + 1

    def __init__(
        self, loc1: Location, loc2: Location, distance: float, id: Optional[int] = None
    ):
        """
        Construct a location edge, which holds attributes that relate two Locations.

        args:
            loc1: the first location
            loc2: the other location
            distance: a measure of distance between the locations
            id: a unique identifier for this edge
        """
        assert loc1 != loc2, "can't have a location self-edge"

        # self.id is unique ID number used to track each edge.
        if id is not None:
            self.id = id
        else:
            self.id = self.next_edge_id

        self.update_id_counter(self.id)

        self.edge = set({loc1, loc2})
        self.distance = distance

        loc1.neighbors.add(loc2.name)
        loc2.neighbors.add(loc1.name)


class Geography:
    def __init__(self, params: ObjMap):
        """
        Umbrella class to initialize/store locations and location edges for a population

        args:
            params: model parameters
        """

        self.locations: Dict[str, Location] = {
            location: Location(location, defn, params)
            for location, defn in params.classes.locations.items()
        }

        self.categories: Dict[str, List[Location]] = {}
        for location in self.locations.values():
            if location.category in self.categories:
                self.categories[location.category].append(location)
            else:
                self.categories[location.category] = [location]

        if params.location.migration.enabled:
            with open(params.location.migration.probs_file, newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    from_loc = row.pop("")
                    prob = float(row.pop("prob", 1))
                    values = list(row.keys())
                    weights = list(map(float, row.values()))
                    assert math.isclose(
                        sum(weights), 1, abs_tol=0.001
                    ), f"Migration weights for {from_loc} must add to 1"
                    if params.location.migration.attribute == "name":
                        self.locations[from_loc].migration_weights["prob"] = prob
                        self.locations[from_loc].migration_weights["weights"] = weights
                        self.locations[from_loc].migration_weights["values"] = values
                    elif params.location.migration.attribute == "category":
                        for location in self.categories[from_loc]:
                            location.migration_weights["prob"] = prob
                            location.migration_weights["weights"] = weights
                            location.migration_weights["values"] = values
                    else:
                        raise ValueError("Unknown migration attribute")

        self.edges: Set[LocationEdge] = set()
        for name, defn in params.location.edges.items():
            if name != "edge_default":
                loc1 = self.locations[defn.location_1]
                loc2 = self.locations[defn.location_2]
                self.edges.add(LocationEdge(loc1, loc2, defn.distance))
