# ============== PARTNERING PROBABILITIES ======================

MSWsexualDurations = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
MSWsexualDurations[1] = {"p_value": (0.27 + 0.22), "min": 1, "max": 6}
MSWsexualDurations[2] = {"p_value": (0.09 + 0.262 + 0.116), "min": 7, "max": 12}
MSWsexualDurations[3] = {"p_value": (0.09 + 0.09), "min": 13, "max": 24}
MSWsexualDurations[4] = {"p_value": (0.09 + 0.09 + 0.07), "min": 25, "max": 36}
MSWsexualDurations[5] = {"min": 37, "max": 48}

# ================ CORE PROBABILITIES ========================


def unsafe_sex(num_acts):
    """
    :Purpose:
        Return probability of unsafe sex in base case
    :Input:
        :num_acts: Number of sex acts
    """
    if num_acts == 0:
        return 0.443
    elif num_acts == 1:
        return 0.481
    elif num_acts < 10:
        return 0.514
    else:
        return 0.759


# TO_REVIEW -these comments?
# F 22% for 6 months or less, 17% for 6mo-1 year, 24% for 1-3 years, 14% for 3-5 years 13% for 5-10years, 10% for 10+ years

# F 74% for 6 months or less, 12% for 6 months to one year, 10% for one to three years, 4% for over three years
HF_jail_duration = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}}
HF_jail_duration[1] = {"p_value": (0.40), "min": 1, "max": 2}
HF_jail_duration[2] = {"p_value": (0.475), "min": 1, "max": 13}
HF_jail_duration[3] = {"p_value": (0.065), "min": 13, "max": 26}
HF_jail_duration[4] = {"p_value": (0.045), "min": 26, "max": 78}
HF_jail_duration[5] = {"p_value": (0.01), "min": 78, "max": 130}
HF_jail_duration[6] = {"p_value": (0.01), "min": 130, "max": 260}

# M 61% for six months or less, 13% for six months to one year, 16% for 1-3 years, 4% for 3-5 years, 4% for 5-10 years

# M 13% for 6 months or less, 8% for 6 mo-1year, 20% for 1-3 years, 11% for 3-5 years, 16% for 5-10 years, 30% for 10+ years
HM_jail_duration = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}}
HM_jail_duration[1] = {"p_value": (0.43), "min": 1, "max": 2}
HM_jail_duration[2] = {"p_value": (0.50), "min": 1, "max": 13}
HM_jail_duration[3] = {"p_value": (0.02), "min": 13, "max": 26}
HM_jail_duration[4] = {"p_value": (0.02), "min": 26, "max": 78}
HM_jail_duration[5] = {"p_value": (0.03), "min": 78, "max": 130}
HM_jail_duration[6] = {"p_value": (0.01), "min": 130, "max": 260}


def adherence_prob(adherenceStat):
    if adherenceStat == 1:
        return 0.0051
    elif adherenceStat == 2:
        return 0.0039
    elif adherenceStat == 3:
        return 0.0032
    elif adherenceStat == 4:
        return 0.0025
    elif adherenceStat == 5:
        return 0.0008
    else:
        return 0.0051


# ======================== HIVABM_Population =========================

# Quick calc based on traj: 0.596 Female -> 0,01175/0.75 HIV -> 0.75 Tested+ = 11.75% HIV DIAG
def get_strat_white(MT_PrEP, MT_Incar, MT_NoIncar):
    """
    Get dictionary of probabilities for white race
    """
    StratW = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}}
    if MT_PrEP:
        StratW["MSM"] = {
            "POP": 1.0,
            "HIV": 0.036195675,
            "AIDS": 0.51,
            "HAARTprev": 0.33,
            "INCARprev": 0.0,
            "TestedPrev": 0.816,
            "mNPart": 3,
        }
        StratW["HM"] = {
            "POP": 0.0,
            "HIV": 0.0,
            "AIDS": 0.0,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.0,
        }
        StratW["HF"] = {
            "POP": 0.0,
            "HIV": 0.0,
            "AIDS": 0.0,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.0,
        }
        StratW["IDU"] = {
            "POP": 0.0,
            "HIV": 0.0,
            "AIDS": 0.0,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.0,
        }
    elif MT_Incar:
        StratW["MSM"] = {
            "POP": 0.0,
            "HIV": 0.05,
            "AIDS": 0.05,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.86,
            "mNPart": 3,
        }
        StratW["HM"] = {
            "POP": 0.4044435412,
            "HIV": 0.0158 / 0.75,
            "AIDS": 0.678,
            "HAARTprev": 0.03748,
            "INCARprev": 0.0279,
            "TestedPrev": 0.75,
            "mNPart": 3,
        }
        StratW["HF"] = {
            "POP": 0.5955564588,
            "HIV": 0.01175 / 0.75,
            "AIDS": 0.573,
            "HAARTprev": 0.04054,
            "INCARprev": 0.0,
            "TestedPrev": 0.75,
            "mNPart": 2,
        }
        StratW["IDU"] = {
            "POP": 0.0,
            "HIV": 0.0,
            "AIDS": 0.0,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.0,
            "mNPart": 3,
        }
    elif MT_NoIncar:
        StratW["MSM"] = {
            "POP": 0.0,
            "HIV": 0.05,
            "AIDS": 0.05,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.86,
            "mNPart": 3,
        }
        StratW["HM"] = {
            "POP": 0.4044435412,
            "HIV": 0.0158 / 0.75,
            "AIDS": 0.678,
            "HAARTprev": 0.03748,
            "INCARprev": 0.0,
            "TestedPrev": 0.75,
            "mNPart": 3,
        }
        StratW["HF"] = {
            "POP": 0.5955564588,
            "HIV": 0.01175 / 0.75,
            "AIDS": 0.573,
            "HAARTprev": 0.04054,
            "INCARprev": 0.0,
            "TestedPrev": 0.75,
            "mNPart": 2,
        }
        StratW["IDU"] = {
            "POP": 0.0,
            "HIV": 0.0,
            "AIDS": 0.0,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.0,
            "mNPart": 3,
        }
    return StratW


def get_strat_black():  # TO_REVIEW - why all zero?
    """
    Return dictionary of probabilities for black race
    """
    StratB = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}}
    StratB["MSM"] = {
        "POP": 0.0,
        "HIV": 0.0,
        "AIDS": 0.0,
        "HAARTprev": 0.0,
        "INCARprev": 0.0,
        "TestedPrev": 0.0,
    }
    StratB["HM"] = {
        "POP": 0.0,
        "HIV": 0.0,
        "AIDS": 0.0,
        "HAARTprev": 0.0,
        "INCARprev": 0.0,
        "TestedPrev": 0.0,
    }
    StratB["HF"] = {
        "POP": 0.0,
        "HIV": 0.0,
        "AIDS": 0.0,
        "HAARTprev": 0.0,
        "INCARprev": 0.0,
        "TestedPrev": 0.0,
    }
    StratB["IDU"] = {
        "POP": 0.0,
        "HIV": 0.0,
        "AIDS": 0.0,
        "HAARTprev": 0.0,
        "INCARprev": 0.0,
        "TestedPrev": 0.0,
    }

    return StratB


def get_IDU():
    """
    Nested dictionary for probability lookup for IDU population
    """
    IDU = {"MSM": {}, "HM": {}, "WSW": {}, "HF": {}}
    IDU["MSM"] = {"HIV": 0.55, "AIDS": 0.058, "HAARTa": 0}
    IDU["HM"] = {"HIV": 0.42, "AIDS": 0.058, "HAARTa": 0}
    IDU["WSW"] = {"HIV": 0.53, "AIDS": 0.058, "HAARTa": 0}
    IDU["HF"] = {"HIV": 0.39, "AIDS": 0.058, "HAARTa": 0}
    return IDU


def get_NIDU():
    """
    Nested dictionary for probability lookup for NIDU population
    """
    NIDU = {"MSM": {}, "HM": {}, "WSW": {}, "HF": {}}
    NIDU["MSM"] = {"HIV": 0.18, "AIDS": 0.02, "HAARTa": 0}
    NIDU["HM"] = {"HIV": 0.048, "AIDS": 0.002, "HAARTa": 0}
    NIDU["WSW"] = {"HIV": 0.048, "AIDS": 0.002, "HAARTa": 0}
    NIDU["HF"] = {"HIV": 0.048, "AIDS": 0.002, "HAARTa": 0}
    return NIDU


def get_ND():
    """
    Nested dictionary for probability lookup for ND population
    """
    ND = {"MSM": {}, "HM": {}, "WSW": {}, "HF": {}}
    ND["Type"] = ([0, 1, 2, 3], [0.469, 0.493, 0.022, 0.016])
    ND["MSM"] = {"HIV": 0.08, "AIDS": 0.02, "HAARTa": 0}
    ND["HM"] = {"HIV": 0.015, "AIDS": 0.0003, "HAARTa": 0}
    ND["WSW"] = {"HIV": 0.012, "AIDS": 0.0003, "HAARTa": 0}
    ND["HF"] = {"HIV": 0.012, "AIDS": 0.0003, "HAARTa": 0}
    return ND


def jail_duration():
    """
    Binned dictionary for jail durations
    """
    jailDuration = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}}
    jailDuration[1] = {"p_value": (0.14), "min": 1, "max": 13}
    jailDuration[2] = {"p_value": (0.09), "min": 13, "max": 26}
    jailDuration[3] = {"p_value": (0.20), "min": 26, "max": 78}
    jailDuration[4] = {"p_value": (0.11), "min": 78, "max": 130}
    jailDuration[5] = {"p_value": (0.16), "min": 130, "max": 260}
    jailDuration[6] = {"p_value": (0.30), "min": 260, "max": 520}
    return jailDuration


def get_mean_num_partners(drug_type, rand_generator):
    diceroll = rand_generator.random()

    if drug_type == "IDU":
        if diceroll < 0.389:
            return 1
        elif diceroll < 0.389 + 0.150:
            return 2
        elif diceroll < 0.389 + 0.150 + 0.089:
            return 3
        elif diceroll < 0.389 + 0.150 + 0.089 + 0.067:
            return 4
        elif diceroll < 0.389 + 0.150 + 0.089 + 0.067 + 0.098:
            return rand_generator.randrange(5, 6, 1)
        elif diceroll < 0.389 + 0.150 + 0.089 + 0.067 + 0.098 + 0.108:
            return rand_generator.randrange(7, 10, 1)
        elif diceroll < 0.389 + 0.150 + 0.089 + 0.067 + 0.098 + 0.108 + 0.02212:
            return rand_generator.randrange(17, 30, 1)
        else:
            return rand_generator.randrange(31, 60, 1)
    else:
        if diceroll < 0.84:
            return 1
        elif diceroll < 0.84 + 0.13:
            return 2
        else:
            return rand_generator.randrange(3, 4, 1)
