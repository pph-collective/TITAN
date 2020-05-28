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
def get_death_rate(hiv, aids, race, haart_adh, death_rate, steps_per_year):
    if hiv:
        if aids:  # AIDS DEATH RATE
            p = death_rate.aids

        elif haart_adh > 1:  # HAART DEATH RATE
            p = death_rate.base

        else:  # HIV+ DEATH RATE
            p = death_rate.hiv

    else:  # NON HIV DEATH RATE
        p = death_rate.base

    # putting it into per 1 person-month from per 1000 person years
    return p / (steps_per_year * 1000.0)
