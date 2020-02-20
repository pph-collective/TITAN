from dotmap import DotMap

params = DotMap()
"""
Settings profiles for different locations
"""


PARTNERTURNOVER = 7.5  # Partner acquisition parameters (higher number more partnering)

cal_NeedleScaling = 1.0  # PWID transmission probability scaling factor
cal_SexualScaling = 1.0  # Sexual transmission probability scaling factor
params.calibration.transmission = (
    1.0  # 0.92 # Global transmission probability scaling factor
)
params.calibration.acute = (
    4.3  # Infectivity multiplier ratio for Acute status infections
)
params.calibration.risk_reduction.transmission = (
    0.25  # Risk reduction in transmission probability for agents diagnosed
)
params.calibration.risk_reduction.haart = (
    3.0  # Scaling factor for effectiveness of ART therapy on xmission P
)
params.calibration.test_frequency = 0.5  # Scaling factor for testing frequency
params.calibration.mortality = 0.05  # Scaling factor for all cause mortality rates
params.calibration.aids_progression = (
    0.05  # Scaling factor for all progression to AIDS from HIV rates
)
params.calibration.haart_coverage = (
    1.45  # Scaling factor for enrollment on ART probability
)
params.calibration.incarceration = 1.0
params.calibration.race_transmission = 2.9

# High risk params
params.high_risk.partner_scale = (
    300  # Linear increase to partner number during HR period
)
params.high_risk.proportion = (
    0.3  # Proportion of people who enter HR group when partner incarcerated
)
params.high_risk.sex_based.HM.duration = 6  # Duration of high risk for males
params.high_risk.sex_based.HF.duration = 6  # Duration of high risk for females


# Misc. params
params.assort_mix.coefficient = (
    0.05  # Proportion of race1 mixing with race2 when partnering.
)
params.needle_exchange.prevalence = 1.0

# Incarceration params
inc_JailMax = 9
inc_JailMin = 1
inc_JailTestProb = 0.69
inc_PrisMax = 60
inc_PrisMin = 6
inc_PrisTestProb = 0.69
inc_PropPrison = 0.5
params.incar.haart.prob = 0.51
params.incar.haart.adherence = 0.21
params.incar.haart.discontinue = 0.12
inc_Recidivism = 0.267

# PrEP params
params.prep.type = "Inj"  # Oral/Inj PrEP modes
params.prep.target = (
    0.05  # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
)
params.prep.start = 0  # Start date for PrEP program (0 for start of model)
params.PWID.prep.adherence = 0.82  # Probability of being adherent
params.prep.efficacy.adherent = 0.96  # Efficacy of adherence PrEP
params.prep.efficacy.non_adherant = 0.76  # Efficacy of non-adherence PrEP
PrEP_falloutT = 0  # During PrEP remains effective post discontinuation
PrEP_resist = 0.0  # 1
params.prep.discontinue = 0.15

if params.prep.type == "Oral":
    params.PWID.prep.adherence = 0.92
    params.prep.efficacy.adherent = 0.96
    params.prep.efficacy.non_adherant = 0.76
    PrEP_falloutT = 1
    params.prep.discontinue = 0.15
elif params.prep.type == "Inj":
    params.PWID.prep.adherence = 1.0
    params.prep.efficacy.adherent = 0.90
    params.prep.efficacy.non_adherant = 1.00
    PrEP_falloutT = 9
    params.prep.discontinue = 0.04


"""
RaceClass is a distinct racial/ethnic/social classification for demographics of the population.
ID is the specific mode of partnership the agent engages in (ie MSM, HM, HF, PWID)
RaceClass agent classifier template
"""
RC_template = {
    "Race": None,  # Race of demographic
    "Class": None,  # Classification of networking
    "POP": 0.0,  # Percentage of total agent population that are ID
    "HIV": 0.0,  # Proportion of total ID population that are HIV
    "AIDS": 0.0,  # Proportion of total HIV_ID that are AIDS
    "HAARTprev": 0.0,  # Proportion of HIV_TESTED_ID that are enrolled on ART
    "INCARprev": 0.0,  # Proportion of ID that are incarcerated
    "TestedPrev": 0.0,  # Proportion of HIV_ID that are tested positively
    "mNPart": 0.0,  # Mean number of sex partners
    "NUMPartn": 0.0,  # Number of partners (redundant)
    "NUMSexActs": 0.0,  # Mean number of sex acts with each partner
    "SAFESEX": 0.0,  # Probability of engaging in safe sex (per act)
    "NEEDLESH": 0.0,  # Probability of sharing syringes during join drug use (per act)
    "HIVTEST": 0.0,  # Probability of testing for HIV
    "INCAR": 0.0,  # Probability of becoming incarcerated (rate)
    "HAARTadh": 0.0,  # Adherence to ART therapy
    "HAARTdisc": 0.0,  # Probability of discontinuing ART therapy
    "PrEPdisc": 0.0,  # Probability of discontinuing PrEP treatment
}

RaceClass1 = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
RaceClass1 = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
RaceClass2 = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
for a in ["MSM", "HM", "HF", "PWID"]:
    RaceClass1[a] = dict(RC_template)
    RaceClass2[a] = dict(RC_template)

RaceClass1["MSM"]["POP"] = 1.0
RaceClass1["MSM"]["HIV"] = 0.4
# StratW['MSM'] = {'POP':0.035, 'HIV':0.132, 'AIDS':0.048, 'HAARTprev':0.57, 'INCARprev':0.005, 'TestedPrev':0.84}

RaceClass1["MSM"] = {
    "POP": 1.00,
    "HIV": 0.132,
    "AIDS": 0.048,
    "HAARTprev": 0.57,
    "INCARprev": 0.005,
    "TestedPrev": 0.84,
    "mNPart": 3,
    "NUMPartn": 1.5,
    "NUMSexActs": 5.0,
    "SAFESEX": 0.43,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.055,
    "INCAR": 0.00,  # 0.00014,
    "HAARTadh": 1.0,  # 0.57,
    "HAARTdisc": 0.008,
    "PrEPdisc": params.prep.discontinue,
}

RaceClass1["PWID"] = {
    "POP": 0.20,  # 0.005,
    "HIV": 0.19,
    "AIDS": 0.075,
    "HAARTprev": 0.0,
    "INCARprev": 0.0,
    "TestedPrev": 0.93,
    "mNPart": 3,
    "NUMPartn": 1.2,
    "NUMSexActs": 5.0,
    "SAFESEX": 0.72,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.037,
    "INCAR": 0.00,  # 0.0050,
    "HAARTadh": 0.43,
    "HAARTdisc": 0.018,
    "PrEPdisc": 0.0,
}
RaceClass1["ALL"] = {"Proportion": 0.571, "HAARTdisc": 0.018, "PrEPdisc": 0.0}

RaceClass2 = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
RaceClass2["MSM"] = {
    "POP": 1.00,  # 0.028,
    "HIV": 0.432,
    "AIDS": 0.048,
    "HAARTprev": 0.47,
    "INCARprev": 0.018,
    "TestedPrev": 0.84,
    "mNPart": 3,
    "NUMPartn": 4.9,
    "NUMSexActs": 6.0,
    "SAFESEX": 0.43,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.06,
    "INCAR": 0.00,  # 0.0011,
    "HAARTadh": 1.0,  # 0.34,
    "HAARTdisc": 0.01,
    "PrEPdisc": params.prep.discontinue,
}

RaceClass2["PWID"] = {
    "POP": 0.20,  # 0.005,
    "HIV": 0.19,
    "AIDS": 0.075,
    "HAARTprev": 0.0,
    "INCARprev": 0.0,
    "TestedPrev": 0.93,
    "mNPart": 3,
    "NUMPartn": 3.5,
    "NUMSexActs": 5.0,
    "SAFESEX": 0.69,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.043,
    "INCAR": 0.00,  # 0.0060,
    "HAARTadh": 0.36,
    "HAARTdisc": 0.023,
    "PrEPdisc": 0.0,
}
RaceClass2["ALL"] = {"Proportion": 0.429, "HAARTdisc": 0.018, "PrEPdisc": 0.0}

params.demographics = {"WHITE": RaceClass1, "BLACK": RaceClass2}


"""
Partnership durations and
"""
params.partnership.sex.duration = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
params.partnership.sex.duration[1] = {"prob": (0.323 + 0.262), "min": 1, "max": 6}
params.partnership.sex.duration[2] = {
    "prob": (0.323 + 0.262 + 0.116),
    "min": 7,
    "max": 12,
}
params.partnership.sex.duration[3] = {
    "prob": (0.323 + 0.262 + 0.116 + 0.121),
    "min": 13,
    "max": 24,
}
params.partnership.sex.duration[4] = {
    "prob": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
params.partnership.sex.duration[5] = {"min": 37, "max": 48}

params.partnership.needle.duration = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
params.partnership.needle.duration[1] = {"prob": 1.0, "min": 1, "max": 6}
params.partnership.needle.duration[2] = {"prob": (0.323 + 0.262), "min": 1, "max": 6}
params.partnership.needle.duration[3] = {"prob": (0.323 + 0.262), "min": 1, "max": 6}
params.partnership.needle.duration[4] = {"prob": (0.323 + 0.262), "min": 1, "max": 6}
params.partnership.needle.duration[5] = {"min": 1, "max": 6}

PartnershipDurations = {
    "SEX": params.partnership.sex.duration,
    "NEEDLE": params.partnership.needle.duration,
}


"""
Sexual and injection transmission probabilities
"""
SexTrans = {"MSM": {}, "HM": {}, "HF": {}}
SexTrans["MSM"] = {
    "0": 0.005,
    "1": 0.005,
    "2": 0.004,
    "3": 0.002,
    "4": 0.001,
    "5": 0.0001,
}
SexTrans["HM"] = {
    "0": 0.001,
    "1": 0.001,
    "2": 0.0008,
    "3": 0.0004,
    "4": 0.0002,
    "5": 0.0001,
}
SexTrans["HF"] = {
    "0": 0.001,
    "1": 0.001,
    "2": 0.0008,
    "3": 0.0004,
    "4": 0.0002,
    "5": 0.0001,
}

params.partnership.needle.transmission = {
    "0": 0.007,
    "1": 0.007,
    "2": 0.0056,
    "3": 0.0028,
    "4": 0.0014,
    "5": 0.0002,
}

TransmissionProbabilities = {
    "SEX": SexTrans,
    "NEEDLE": params.partnership.needle.transmission,
}
