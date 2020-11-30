from . import utils
from .location import Location

# ================ CORE PROBABILITIES ========================


@utils.memo
def adherent_prob(adherent: bool) -> float:
    """
    Mapping from HAART adherence levels to probabilities.

    args:
        adherent: Whether an agent is HAART adherent

    returns:
        probability of agent transitioning from HIV+ to AIDS
    """
    if adherent:
        return 0.0008
    else:
        return 0.00368


@utils.memo
def get_death_rate(
    hiv: bool,
    aids: bool,
    drug_type: str,
    haart_adh: int,
    race: str,
    location: Location,
    steps_per_year: int,
) -> float:
    """
    Find the death rate of an agent given a set of attributes.

    args:
        hiv: whether the agent is HIV+
        aids: whether the agent has AIDS
        drug_type: whether the PWID base death rate should be used or the base one
        haart_adh: whether an agent is haart adherent
        race: the race of the agent
        location: agent's location
        steps_per_year: the number of model steps in a year

    returns:
        the probability of an agent with these characteristics dying in a given time step
    """
    param = location.params.demographics

    if drug_type == "Inj":
        death_param = param[race].PWID.death_rate
    else:
        death_param = param[race].death_rate

    p = death_param.base

    if aids:
        p *= death_param.aids
    elif hiv:
        if haart_adh:
            p *= death_param.haart_adherent
        else:
            p *= death_param.hiv

    # putting it into per 1 person-month from per 1000 person years
    return p / (steps_per_year * 1000.0)
