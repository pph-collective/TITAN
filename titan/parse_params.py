import paraml  # type: ignore

import os
from inspect import getsourcefile
from pathlib import Path
import shutil
import math
from typing import Optional, Dict
from copy import deepcopy


class ObjMap(dict):
    """
    A dictionary-like class which allows accessing members either using standard
    dictionary notation or dots.  Note the hash function is hard-coded - beware.
    """

    def __init__(self, d: Dict):
        for k, v in d.items():
            if isinstance(v, dict):
                v = self.__class__(v)
            self[k] = v

    def __getattribute__(self, k):
        try:
            return self[k]
        except KeyError:
            return object.__getattribute__(self, k)

    def __setattr__(self, k, v):
        return self.__setitem__(k, v)

    def __hash__(self):
        return 1234567890

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, state):
        self.__dict__.update(state)

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.items():
            result[k] = deepcopy(v, memo)
        return result


# ============== PARSING FUNCTIONS ======================


def check_params(params: ObjMap):
    """
    Consistency checks for param populations
    """
    race_pop = 0
    for race in params.classes.races:
        r_dems = params.demographics[race]
        race_pop += r_dems.ppl
        sex_type_pop = 0
        for st, st_dems in r_dems.sex_type.items():
            if st in list(params.classes.sex_types.keys()):
                sex_type_pop += st_dems.ppl

            drug_type_pop = 0
            for dt, dt_dems in st_dems.drug_type.items():
                if dt in list(params.classes.drug_types):
                    drug_type_pop += dt_dems.ppl

            assert math.isclose(
                drug_type_pop, 1, abs_tol=0.001
            ), f"ppl of {race}'s {st}'s drug_types must add to 1. Currently adding to {drug_type_pop}"

        assert math.isclose(
            sex_type_pop, 1, abs_tol=0.001
        ), f"ppl of {race}'s sex_types must add to 1. Currently adding to {sex_type_pop}"

    assert math.isclose(race_pop, 1, abs_tol=0.001), "ppl of races must add to 1"

    loc_pop = 0
    for location in params.classes.locations.values():
        loc_pop += location.ppl

    assert math.isclose(loc_pop, 1, abs_tol=0.001), "ppl of locations must add to 1"

    for param, assort in params.assort_mix.items():
        assort_value = 0
        for ptnr_value in assort.partner_values.values():
            assort_value += ptnr_value
        assert math.isclose(
            assort_value, 1, abs_tol=0.001
        ), f"assort values must add to 1, not {assort_value} in {param}"


def create_params(
    setting_name: Optional[str],
    param_path: str,
    outdir: str,
    error_on_unused: bool = False,
) -> ObjMap:
    """
    Entry function - given the path to the setting, params, output directory and whether
    or not to use the base setting. Parse and create a params (ObjMap) object.

    args:
        setting_name: path to a settings file or directory or `None`
        param_path: path to parameter file or directory
        outdir: path to directory where computed params will be saved
        error_on_unused: throw a hard error if there are unused parameters, otherwise warnings are only printed

    returns:
        computed/validated model paramters with defaults filled in where needed
    """
    # find defs, where we are in the code for settings and base
    filename = getsourcefile(create_params)  # what is the sourcefile for this function
    if filename is not None:
        parent = Path(filename).resolve().parent
    else:
        raise Exception("can't find where I am in the code?")

    param_defs = os.path.join(parent, "params")

    param_paths = []

    # merge setting and params
    if setting_name is not None:
        # check if it's a known setting or pass it through as a path
        if setting_name in os.listdir(os.path.join(parent, "settings")):
            param_paths.append(os.path.join(parent, "settings", setting_name))
        else:
            param_paths.append(setting_name)

    param_paths.append(param_path)

    parsed = paraml.create_params(
        param_defs,
        *param_paths,
        out_path=os.path.join(outdir, "params.yml"),
        error_on_unused=error_on_unused,
    )

    parsed = ObjMap(parsed)
    check_params(parsed)

    # copy migration file if enabled
    if parsed.location.migration.enabled:
        shutil.copy(
            parsed.location.migration.probs_file,
            os.path.join(outdir, "migration_probs.csv"),
        )

    return parsed
