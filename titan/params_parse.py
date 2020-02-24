import oyaml as yaml  # type: ignore
import os
import sys
import collections
from inspect import currentframe, getframeinfo
from pathlib import Path

from dotmap import DotMap  # type: ignore


def check_item(val, d, keys=None):
    if "min" in d:
        assert val >= d["min"]
    if "max" in d:
        assert val <= d["max"]
    if d["type"] == "int":
        isinstance(val, int)
    if d["type"] == "float":
        if isinstance(val, int):
            val = float(val)
        assert isinstance(val, float)
    if d["type"] == "boolean":
        assert isinstance(val, bool)
    if d["type"] == "enum":
        assert val in d["values"]
    if d["type"] == "array":
        assert all(x in d["values"] for x in val)
    if d["type"] == "keys":
        assert all(x in keys for x in val)
    return val


def get_item(key, d, param):
    if key in param:
        val = param[key]
        return check_item(val, d)
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

    parsed_bins = {}
    for bin, val in bins.items():
        try:
            int(bin)
        except:
            print("Bins must be integers")
            raise

        for field, defn in d["fields"].items():
            assert field in val
            val[field] = check_item(val[field], defn)

        parsed_bins[int(bin)] = val

    return parsed_bins


def get_defn(key, d, param):
    if key not in param:
        parsed = d["default"]
    else:
        parsed = param[key]

    # check definitions
    for k, val in parsed.items():
        for field, defn in d["fields"].items():
            assert field in val
            val[field] = check_item(val[field], defn, parsed.keys())

    return parsed


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
            elif v["type"] == "definition":
                parsed[k] = get_defn(k, v, params)
            else:
                parsed[k] = get_item(k, v, params)
        else:
            parsed[k] = parse_params(v, params.get(k, {}), pops)

    return parsed


def parse_classes(defs, params):
    # add sex types to populations
    if "sex_types" in params.get("classes", []):
        params["classes"]["populations"] = (
            params["classes"].get(
                "populations", defs["classes"]["populations"]["default"]
            )
            + list(params["classes"]["sex_types"].keys())
        )

    sex_type_keys = list(defs["classes"]["sex_types"]["default"].keys())

    defs["classes"]["populations"]["default"] += sex_type_keys
    defs["classes"]["populations"]["values"] += sex_type_keys

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


def build_yaml(path):
    yml = {}
    if os.path.isdir(path):
        for file in os.listdir(path):
            with open(os.path.join(path, file)) as f:
                this_yml = yaml.safe_load(f)
                yml.update(this_yml)
    else:
        with open(path) as f:
            this_yml = yaml.safe_load(f)
            yml.update(this_yml)

    return yml


def create_params(setting_path, param_path, outdir):
    filename = getframeinfo(currentframe()).filename
    parent = Path(filename).resolve().parent
    root = os.path.join(parent, "params")

    defs = build_yaml(root)
    params = build_yaml(param_path)

    # merge setting and params
    if setting_path is not None:
        setting = build_yaml(setting_path)
        params = merge(setting, params)

    pops = parse_classes(defs, params)
    parsed = parse_params(defs, params, pops)

    with open(os.path.join(outdir, "params.yml"), "w") as f:
        yaml.dump(parsed, f)

    # this is just for dev why params change - TODO: DELETE
    with open("dotmap_params.txt", "w") as f:
        f.write("dot,old\n")
        print_dotmap(parsed, "", f)

    return DotMap(parsed)


if __name__ == "__main__":
    create_params(None, "tests/params/basic.yml", "results")
