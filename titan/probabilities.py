from . import utils

# ================ CORE PROBABILITIES ========================


@utils.memo
def adherence_prob(adherence):
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
def get_death_rate(hiv, aids, drug_use, haart_adh, param, steps_per_year):
    if drug_use == "Inj":
        death_param = param.PWID.death_rate
    else:
        death_param = param.death_rate
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
