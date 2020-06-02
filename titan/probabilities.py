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
def get_death_rate(hiv, aids, so, drug_use, haart_adh, param, steps_per_year):
    if drug_use == "Inj":
        p = param.PWID.death_rate.base
        if aids:
            p *= param.PWID.death_rate.aids
        elif hiv:
            if haart_adh == 5:
                p *= param.PWID.death_rate.haart_adherent
            else:
                p *= param.PWID.death_rate.hiv
    else:
        p = param.death_rate.base
        if hiv:
            if aids:  # AIDS DEATH RATE
                p *= param.death_rate.aids

            elif haart_adh == 5:  # HAART DEATH RATE
                p *= param.death_rate.haart_adherent

            else:  # HIV+ DEATH RATE
                p = param.death_rate.hiv

    # putting it into per 1 person-month from per 1000 person years
    return p / (steps_per_year * 1000.0)
