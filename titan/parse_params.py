import oyaml as yaml  # type: ignore
import os
import collections
from inspect import currentframe, getframeinfo
from pathlib import Path
import math


class ObjMap(dict):
    """
    A dictionary-like class which allows accessing members either using standard
    dictionary notation or dots.  Note the hash function is hard-coded - beware.
    """

    def __init__(self, d):
        for k, v in d.items():
            if isinstance(v, dict):
                v = self.__class__(v)
            self[k] = v

    def __getattr__(self, k):
        return self.__getitem__(k)

    def __hash__(self):
        return 1234567890

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__.update(state)


# ============== PARSING FUNCTIONS ======================


def check_item(val, d, keys=None, pops=None):
    """
    Checks if an item meets the requirements of the field's definition.
    """
    if "min" in d:
        assert val >= d["min"], f"{val} must be greater than {d['min']}"
    if "max" in d:
        assert val <= d["max"], f"{val} must be less than {d['max']}"
    if d["type"] == "int":
        assert isinstance(val, int), f"{val} must be an integer"
    if d["type"] == "float":
        if isinstance(val, int):
            val = float(val)
        assert isinstance(val, float), f"{val} must be a float"
    if d["type"] == "boolean":
        assert isinstance(val, bool), f"{val} must be a bool"
    if d["type"] == "enum":
        if "values" in d:
            values = d["values"]
        elif "class" in d:
            values = pops[d["class"]]
        assert val in values, f"{val} not in {values}"
    if d["type"] == "array":
        if "values" in d:
            values = d["values"]
        elif "class" in d:
            values = pops[d["class"]]
        assert isinstance(val, list), f"{val} must be an array"
        assert all(x in values for x in val), f"{val} not in {values}"
    if d["type"] == "keys":
        assert isinstance(val, list), f"{val} must be an array of keys"
        assert all(x in keys for x in val), f"{keys} not in {keys}"
    return val


def get_item(key, d, param, pops=None):
    """
    Get and check item from the params, falling back on the definitions default.
    """
    if key in param:
        val = param[key]
        return check_item(val, d, pops=pops)
    else:
        return d["default"]


def merge(d1, d2):
    """return new merged dict of dicts"""
    if isinstance(d1, collections.abc.Mapping) and isinstance(
        d2, collections.abc.Mapping
    ):
        for k, v in d1.items():
            if k in d2:
                d2[k] = merge(v, d2[k])
        d3 = d1.copy()
        d3.update(d2)
        return d3
    else:
        return d2


def get_bins(key, d, param, pops):
    """
    Get and validate a type == bin definition
    """
    if key not in param:
        return d["default"]

    bins = merge(d["default"], param[key])

    parsed_bins = {}
    for bin, val in bins.items():
        try:
            int(bin)
        except ValueError:
            print("Bins must be integers")
            raise

        for field, defn in d["fields"].items():
            assert field in val, f"{field} must be in {val}"
            val[field] = check_item(val[field], defn, pops=pops)

        parsed_bins[int(bin)] = val

    return parsed_bins


def get_defn(key, d, param, pops):
    """
    Get and validate a type == definition definition
    """
    if key not in param:
        parsed = d["default"]
    else:
        parsed = param[key]

    # check definitions
    for k, val in parsed.items():
        for field, defn in d["fields"].items():
            assert field in val, f"{field} must be in {val}"
            val[field] = check_item(val[field], defn, keys=parsed.keys(), pops=pops)

    return parsed


def parse_params(defs, params, pops):
    """
    Recursively parse the passed params, using the definitions to validate
    and provide defaults.
    """
    parsed = {}
    # params is a scalar, return it
    if not isinstance(params, dict):
        return params

    # handles case of bin as direct default item
    if "default" in defs and defs["type"] == "bin":
        return get_bins("dummy", defs, {"dummy": params}, pops)

    for k, v in defs.items():
        # assumes all v are dicts, as otherwise it would have returned
        if "default" in v:
            if v["type"] == "sub-dict":
                parsed[k] = {}
                field = v["keys"].pop(0)
                for val in pops[field]:
                    parsed[k][val] = parse_params(
                        v["default"], params.get(k, {}).get(val, {}), pops
                    )

                    if len(v["keys"]) > 0:
                        field2 = v["keys"][0]
                        for val2 in pops[field2]:
                            parsed[k][val][val2] = parse_params(
                                v["default"],
                                params.get(k, {}).get(val, {}).get(val2, {}),
                                pops,
                            )

            elif v["type"] == "bin":
                parsed[k] = get_bins(k, v, params, pops)
            elif v["type"] == "definition":
                parsed[k] = get_defn(k, v, params, pops)
            else:
                parsed[k] = get_item(k, v, params, pops)
        else:
            parsed[k] = parse_params(v, params.get(k, {}), pops)

    return parsed


def parse_classes(defs, params):
    """
    Parse the classes definition first as it is needed in parsing the full params.
    """
    # add sex types to populations
    if "sex_types" in params.get("classes", []):
        params["classes"]["populations"] = params["classes"].get(
            "populations", defs["classes"]["populations"]["default"]
        ) + list(params["classes"]["sex_types"].keys())

    sex_type_keys = list(defs["classes"]["sex_types"]["default"].keys())

    defs["classes"]["populations"]["default"] += sex_type_keys
    defs["classes"]["populations"]["values"] += sex_type_keys

    return parse_params(defs["classes"], params.get("classes", {}), {})


def build_yaml(path):
    """
    Read in a yaml or folder of yamls into a dictionary.
    """
    yml = {}
    if os.path.isdir(path):
        for file in os.listdir(path):
            if ".yml" in file or ".yaml" in file:
                with open(os.path.join(path, file)) as f:
                    this_yml = yaml.safe_load(f)
                    yml.update(this_yml)
    else:
        with open(path) as f:
            assert ".yml" in path or ".yaml" in path
            this_yml = yaml.safe_load(f)
            yml.update(this_yml)

    return yml


def check_params(params):
    """
    Consistency checks for param populations
    """
    race_pop = 0

    for race in params.classes.races:
        r_dems = params.demographics[race]
        race_pop += r_dems.ppl
        sex_type_pop = 0
        for st, st_dems in r_dems.items():
            if st in list(params.classes.sex_types.keys()):
                sex_type_pop += st_dems.ppl

        assert math.isclose(
            sex_type_pop, 1, abs_tol=0.001
        ), f"ppl of {race}'s sex_types must add to 1"

    assert math.isclose(race_pop, 1, abs_tol=0.001), f"ppl of races must add to 1"


def create_params(setting_path, param_path, outdir, use_base=True):
    """
    Entry funtion - given the path to the setting, params, output directory and
    whether or not to use the base setting. Parse and create a params (ObjMap) object.
    """
    filename = getframeinfo(currentframe()).filename
    parent = Path(filename).resolve().parent
    root = os.path.join(parent, "params")

    defs = build_yaml(root)
    params = build_yaml(param_path)

    # merge setting and params
    if setting_path is not None:
        setting = build_yaml(setting_path)
        params = merge(setting, params)

    if use_base:
        base_dir = os.path.join(parent, "..", "settings", "base")
        base = build_yaml(base_dir)
        params = merge(base, params)

    pops = parse_classes(defs, params)
    parsed = parse_params(defs, params, pops)

    with open(os.path.join(outdir, "params.yml"), "w") as f:
        yaml.dump(parsed, f)

    parsed = ObjMap(parsed)
    check_params(parsed)

    return parsed
