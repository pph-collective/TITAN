from . import utils
from .location import Location

# ================ CORE PROBABILITIES ========================


@utils.memo
def get_death_rate(
    hiv: bool,
    aids: bool,
    drug_type: str,
    sex_type: str,
    haart_adh: bool,
    race: str,
    location: Location,
    steps_per_year: int,
    exit_class: str,
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
        exit_class: the exit class to access in agent params

    returns:
        the probability of an agent with these characteristics dying in a given time step
    """
    param = location.params.demographics

    death_param = param[race].sex_type[sex_type].drug_type[drug_type].exit[exit_class]

    p = death_param.base

    if aids:
        p *= death_param.aids
    elif hiv:
        if haart_adh:
            p *= death_param.haart_adherent
        else:
            p *= death_param.hiv

    # putting it into per 1 person-month from per 1000 person years
    return p / (1000 * steps_per_year)
