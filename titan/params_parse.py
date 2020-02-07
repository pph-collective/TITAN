import oyaml as yaml  # type: ignore
import os
import sys
import collections
from inspect import currentframe, getframeinfo
from pathlib import Path

from dotmap import DotMap


def get_item(key, d, param):
    if key in param:
        val = param[key]
        if "min" in d:
            assert val >= d["min"]
        if "max" in d:
            assert val <= d["max"]
        if d["type"] == "enum":
            assert val in d["values"]
        if d["type"] == "array":
            assert all(x in d["values"] for x in val)
        return val
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


def get_bins(key, d, param):
    if key not in param:
        return d["default"]

    bins = merge(d["default"], param[key])

    for bin, val in bins.items():
        try:
            int(bin)
        except:
            print("Bins must be integers")
            raise

        for field, defn in d["fields"].items():
            assert field in val
            if "min" in defn:
                assert val[field] >= defn["min"]

            if "max" in defn:
                assert val[field] <= defn["max"]

    return bins


def parse_params(defs, params, pops):
    parsed = {}

    # params is a scalar, return it
    if not isinstance(params, dict):
        return params

    # handles case of bin as direct default item
    if "default" in defs and defs["type"] == "bin":
        return get_bins("dummy", defs, {"dummy": params})

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
                parsed[k] = get_bins(k, v, params)
            else:
                parsed[k] = get_item(k, v, params)
        else:
            parsed[k] = parse_params(v, params.get(k, {}), pops)

    return parsed


def parse_classes(defs, params):
    # add sex types to popultaions
    if "sex_types" in params["classes"]:
        params["classes"]["populations"] = (
            params["classes"].get(
                "populations", defs["classes"]["populations"]["default"]
            )
            + params["classes"]["sex_types"]
        )

    defs["classes"]["populations"]["default"] += defs["classes"]["sex_types"]["default"]
    defs["classes"]["populations"]["values"] += defs["classes"]["sex_types"]["values"]

    return parse_params(defs["classes"], params.get("classes", {}), {})

def print_dotmap(params, prefix, file_handle):
    for k, v in params.items():
        if prefix == "":
            new_prefix = k
        else:
            new_prefix = f"{prefix}.{k}"
        file_handle.write(f"{new_prefix},\n")
        if isinstance(v, dict):
            print_dotmap(v, new_prefix, file_handle)


def create_params(setting_dict, param_path):
    defs = {}
    filename = getframeinfo(currentframe()).filename
    parent = Path(filename).resolve().parent
    root = os.path.join(parent, "params")
    for file in os.listdir(root):
        with open(os.path.join(root, file), "r") as f:
            this_defs = yaml.safe_load(f)
            defs.update(this_defs)

    with open(param_path, "r") as f:
        params = yaml.safe_load(f)

    # merge setting and params
    params = merge(setting_dict, params)

    pops = parse_classes(defs, params)
    parsed = parse_params(defs, params, pops)

    with open("results/params.yml", "w") as f:
        yaml.dump(parsed, f)

    with open("dotmap_params.txt", "w") as f:
        f.write("dot,old\n")
        print_dotmap(parsed, "", f)

    return DotMap(parsed)


if __name__ == "__main__":
    create_params({}, "tests/params/basic.yml")
