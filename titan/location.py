from typing import Optional, Set
from copy import deepcopy

from .parse_params import ObjMap
from . import utils


class Location:
    """
    :Purpose:
        This class constructs and represents a location within the model.  A location
        can have an arbitrary geographic granularity.

    :Input:
        location_def : ObjMap - Definition of location from the params
    """

    next_location_id = 0

    @classmethod
    def update_id_counter(cls, last_id: int):
        cls.next_location_id = last_id + 1

    def __init__(
        self, name: str, defn: ObjMap, params: ObjMap, id: Optional[int] = None
    ):
        """
        Initialize location object
        """

        # self.id is unique ID number used to track each location.
        if id is not None:
            self.id = id
        else:
            self.id = self.next_location_id

        self.update_id_counter(self.id)

        # location properties
        self.name = name
        self.params = self.scale_params(params)
        self.ppl = defn.ppl  # percent of overall population assigned to this location

        # value/weight maps needed for creating new agents in this location
        self.pop_weights: Dict[str, Dict[str, List[Any]]] = {}
        self.role_weights: Dict[str, Dict] = {}
        self.drug_weights: Dict[str, Dict] = {}
        self.init_weights()

        self.neighbors = set({})  # or maybe edges instead

    def __str__(self):
        return self.name

    def scale_params(self, params: ObjMap):
        """
        Scale the generic parameters with any location based scaling from params.location.scaling
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
        for race in self.params.classes.races:
            self.role_weights[race] = {}
            self.drug_weights[race] = {}
            self.pop_weights[race] = {}
            self.pop_weights[race]["values"] = []
            self.pop_weights[race]["weights"] = []
            for st in self.params.classes.sex_types:
                self.pop_weights[race]["values"].append(st)
                self.pop_weights[race]["weights"].append(
                    self.params.demographics[race][st].ppl
                )
                self.role_weights[race][st] = {}
                self.role_weights[race][st]["values"] = []
                self.role_weights[race][st]["weights"] = []
                self.drug_weights[race][st] = {}
                self.drug_weights[race][st]["values"] = []
                self.drug_weights[race][st]["weights"] = []
                for role, prob in self.params.demographics[race][
                    st
                ].sex_role.init.items():
                    self.role_weights[race][st]["values"].append(role)
                    self.role_weights[race][st]["weights"].append(prob)
                for use_type, prob in self.params.demographics[race][
                    st
                ].drug_type.items():
                    self.drug_weights[race][st]["values"].append(use_type)
                    self.drug_weights[race][st]["weights"].append(prob.init)


class LocationEdge:  # is this a directed or undirected edge?

    @classmethod
    def update_id_counter(cls, last_id: int):
        cls.next_edge_id = last_id + 1

    def __init__(self, loc1: Location, loc2: Location, distance: float, id: Optional[int] = None):
        assert loc1 != loc2, "can't have a location self-edge"

        # self.id is unique ID number used to track each edge.
        if id is not None:
            self.id = id
        else:
            self.id = self.next_edge_id

        self.update_id_counter(self.id)

        self.edge = set({loc1, loc2})
        self.distance = distance


class Geography:
    """
    Umbrella class to initialize/store locations and location edges
    """

    def __init__(self, params: ObjMap):

        self.locations: Dict[str, Location] = {
            location: Location(location, defn, params)
            for location, defn in params.classes.locations.items()
        }

        self.edges: Set[LocationEdge] = set()
        for name, defn in params.location.edges.items():
            if name != "edge_default":
                loc1 = self.locations[defn.location_1]
                loc2 = self.locations[defn.location_2]
                self.edges.add(LocationEdge(loc1, loc2, defn.distance))
