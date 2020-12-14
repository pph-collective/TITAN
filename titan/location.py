from typing import Optional, Set, Dict, List, Any
from copy import deepcopy

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

        # value/weight maps needed for creating new agents in this location
        self.pop_weights: Dict[str, Dict[str, List[Any]]] = {}
        self.role_weights: Dict[str, Dict] = {}
        self.drug_weights: Dict[str, Dict] = {}
        self.init_weights()

        self.edges: Set["LocationEdge"] = set({})  # or maybe edges instead

    def __str__(self):
        return self.name

    def __repr__(self):
        return f"'{self.name}'"

    def __eq__(self, other):
        return self.name == other.name

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
        for race in self.params.classes.races:
            self.role_weights[race] = {}
            self.drug_weights[race] = {}
            self.pop_weights[race] = {}
            self.pop_weights[race]["values"] = []
            self.pop_weights[race]["weights"] = []
            race_param = self.params.demographics[race]
            for st in self.params.classes.sex_types:
                st_param = race_param.sex_type[st]
                self.pop_weights[race]["values"].append(st)
                self.pop_weights[race]["weights"].append(st_param.ppl)
                self.role_weights[race][st] = {}
                self.role_weights[race][st]["values"] = []
                self.role_weights[race][st]["weights"] = []
                self.drug_weights[race][st] = {}
                self.drug_weights[race][st]["values"] = []
                self.drug_weights[race][st]["weights"] = []
                for role, prob in st_param.sex_role.init.items():
                    self.role_weights[race][st]["values"].append(role)
                    self.role_weights[race][st]["weights"].append(prob)
                for dt in self.params.classes.drug_types:
                    dt_param = st_param.drug_type[dt]
                    self.drug_weights[race][st]["values"].append(dt)
                    self.drug_weights[race][st]["weights"].append(dt_param.ppl)


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

        self.edges: Set[LocationEdge] = set()
        for name, defn in params.location.edges.items():
            if name != "edge_default":
                loc1 = self.locations[defn.location_1]
                loc2 = self.locations[defn.location_2]
                self.edges.add(LocationEdge(loc1, loc2, defn.distance))
