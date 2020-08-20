from . import utils
from .parse_params import ObjMap

# ================ CORE PROBABILITIES ========================


@utils.memo
def adherence_prob(adherence: int):
    """
    Mapping from HAART adherence levels to probabilities.

    args:
        adherence: HAART adherence level

    returns:
        probability of agent transitioning from HIV+ to AIDS
    """
    if adherence == 1:
        return 0.0051
    elif adherence == 2:
        return 0.0039
    elif adherence == 3:
        return 0.0032
    elif adherence == 4:
        return 0.0025
    elif adherence == 5:
        return 0.0008
    else:
        return 0.0051


@utils.memo
def get_death_rate(
    hiv: bool,
    aids: bool,
    drug_type: str,
    haart_adh: int,
    race: str,
    param: ObjMap,
    steps_per_year: int,
):
    """
    Find the death rate of an agent given a set of attributes.

    args:
        hiv: whether the agent is HIV+
        aids: whether the agent has AIDS
        drug_type: whether the PWID base death rate should be used or the base one
        haart_adh: level of HAART adherence
        race: the race of the agent
        param: demographic params
        steps_per_year: the number of model steps in a year
    """
    if drug_type == "Inj":
        death_param = param[race].PWID.death_rate
    else:
        death_param = param[race].death_rate

    p = death_param.base

    if aids:
        p *= death_param.aids
    elif hiv:
        if haart_adh == 5:
            p *= death_param.haart_adherent
        else:
            p *= death_param.hiv

    # putting it into per 1 person-month from per 1000 person years
    return p / (steps_per_year * 1000.0)
