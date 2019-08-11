__author__ = "MaximilianKing"


"""
Global variables
    flag_incar      Incarceration effects
    flag_PrEP       PrEP enrollment
    flag_HR         High risk behavior for incar or genPop
    flag_ART        ART therapy enrollment
    flag_DandR      Die and replace functionality

"""

####################
PROCESSES = 1  # number of processes in parallel (quadcore)
N_MC = 1  # total number of iterations (Monte Carlo runs)
N_POP = 23815  # population size
TIME_RANGE = 1  # total time steps to iterate

PrEPTarget = 1.2  # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
PrEPstartT = 48
model = "Incar"  # Model Type for fast flag toggling
####################


"""
Calibration scaling parameters for fitting to empirical data

"""

PARTNERTURNOVER = 7.5  # Partner acquisition parameters (higher number more partnering)

cal_NeedleScaling = 1.0  # IDU transmission probability scaling factor
cal_SexualScaling = 0.8  # Sexual transmission probability scaling factor
cal_pXmissionScaling = 0.92  # Global transmission probability scaling factor


cal_AcuteScaling = 1.0  # Infectivity multiplier ratio for Acute status infections
cal_RR_Dx = 0.5  # Risk reduction in transmission probability for agents diagnosed
cal_RR_HAART = 1.0  # Scaling factor for effectiveness of ART therapy on xmission P
cal_TestFreq = 0.5  # Scaling factor for testing frequency
cal_Mortality = 5.0  # Scaling factor for all cause mortality rates
cal_PROGAIDS = 1.0  # Scaling factor for all progression to AIDS from HIV rates
cal_ART_cov = 0.25  # Scaling factor for enrollment on ART probability

HR_partnerScale = 300  # Linear increase to partner number during HR period
HR_proportion = 0.3  # Proportion of people who enter HR group when partner incarcerated
HR_M_dur = 6  # Duration of high risk for males
HR_F_dur = 6  # Duration of high risk for females


"""
Model Type for fast flag toggling
    flag_incar      Incarceration effects
    flag_PrEP       PrEP enrollment
    flag_HR         High risk behavior for incar or genPop
    flag_ART        ART therapy enrollment
    flag_DandR      Die and replace functionality

"""

####################

####################

if model == "PrEP":
    flag_incar = False
    flag_PrEP = True
    flag_HR = False
    flag_ART = True
    flag_DandR = True
elif model == "Incar":
    flag_incar = True
    flag_PrEP = False
    flag_HR = True
    flag_ART = True
    flag_DandR = True
elif model == "NoIncar":
    flag_incar = False
    flag_PrEP = False
    flag_HR = True
    flag_ART = True
    flag_DandR = True


"""
RaceClass is a distinct racial/ethnic/social classification for demographics of the population.
ID is the specific mode of partnership the agent engages in (ie MSM, HM, HF, PWID)

RaceClass#['ID'] = {'POP':          Percentage of total agent population that are ID
                    'HIV':          Proportion of total ID population that are HIV
                    'AIDS':0.0,     Proportion of total HIV_ID that are AIDS
                    'HAARTprev':    Proportion of HIV_TESTED_ID that are enrolled on ART
                    'INCARprev':    Proportion of ID that are incarcerated
                    'TestedPrev':   Proportion of HIV_ID that are tested positively
                    'mNPart':       Mean number of sex partners
                    'NUMPartn':     Number of partners (redundant)
                    'NUMSexActs':   Mean number of sex acts with each partner
                    'UNSAFESEX':    Probability of engaging in unsafe sex (per act)
                    'NEEDLESH':     Probability of sharing syringes during join drug use (per act)
                    'HIVTEST':      Probability of testing for HIV
                    'INCAR':        Probability of becoming incarcerated (rate)
                    'HAARTadh':     Adherence to ART therapy
                    'HAARTdisc':    Probability of discontinuing ART therapy
                    }
"""


RaceClass1 = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
RaceClass1["MSM"] = {
    "POP": 1.0,
    "HIV": 0.036195675,
    "AIDS": 0.51,
    "HAARTprev": 0.33,
    "INCARprev": 0.0,
    "TestedPrev": 0.816,
    "mNPart": 3,
    "NUMPartn": 1.5,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.43,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.055,
    "INCAR": 0.0,
    "HAARTadh": 0.82,
    "HAARTdisc": 0.08,
}
RaceClass1["HM"] = {
    "POP": 0.0,
    "HIV": 0.0,
    "AIDS": 0.0,
    "HAARTprev": 0.0,
    "INCARprev": 0.0,
    "TestedPrev": 0.0,
    "mNPart": 3,
    "NUMPartn": 1.2,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.83,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.013,
    "INCAR": 0.0,
    "HAARTadh": 0.56,
    "HAARTdisc": 0.01,
}
RaceClass1["HF"] = {
    "POP": 0.0,
    "HIV": 0.0,
    "AIDS": 0.0,
    "HAARTprev": 0.0,
    "INCARprev": 0.0,
    "TestedPrev": 0.0,
    "mNPart": 3,
    "NUMPartn": 1.2,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.83,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.015,
    "INCAR": 0.0000,
    "HAARTadh": 0.53,
    "HAARTdisc": 0.013,
}
RaceClass1["PWID"] = {
    "POP": 0.0,
    "HIV": 0.0,
    "AIDS": 0.0,
    "HAARTprev": 0.0,
    "INCARprev": 0.0,
    "TestedPrev": 0.0,
    "mNPart": 3,
    "NUMPartn": 1.2,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.72,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.037,
    "INCAR": 0.0050,
    "HAARTadh": 0.52,
    "HAARTdisc": 0.018,
}
RaceClass1["ALL"] = {"Proportion": 1.0, "HAARTdisc": 0.018}

RaceClass2 = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
RaceClass2["MSM"] = {
    "POP": 0.0,
    "HIV": 0.0,
    "AIDS": 0.0,
    "HAARTprev": 0.0,
    "INCARprev": 0.0,
    "TestedPrev": 0.0,
    "mNPart": 3,
    "NUMPartn": 4.9,
    "NUMSexActs": 6.0,
    "UNSAFESEX": 0.43,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.06,
    "INCAR": 0.0011,
    "HAARTadh": 0.41,
    "HAARTdisc": 0.01,
}
RaceClass2["HM"] = {
    "POP": 0.0,
    "HIV": 0.0,
    "AIDS": 0.0,
    "HAARTprev": 0.0,
    "INCARprev": 0.0,
    "TestedPrev": 0.0,
    "mNPart": 3,
    "NUMPartn": 3.5,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.70,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.023,
    "INCAR": 0.0011,
    "HAARTadh": 0.34,
    "HAARTdisc": 0.013,
}
RaceClass2["HF"] = {
    "POP": 0.0,
    "HIV": 0.0,
    "AIDS": 0.0,
    "HAARTprev": 0.0,
    "INCARprev": 0.0,
    "TestedPrev": 0.0,
    "mNPart": 3,
    "NUMPartn": 3.5,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.70,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.033,
    "INCAR": 0.00016,
    "HAARTadh": 0.50,
    "HAARTdisc": 0.017,
}
RaceClass2["PWID"] = {
    "POP": 0.0,
    "HIV": 0.0,
    "AIDS": 0.0,
    "HAARTprev": 0.0,
    "INCARprev": 0.0,
    "TestedPrev": 0.0,
    "mNPart": 3,
    "NUMPartn": 3.5,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.69,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.043,
    "INCAR": 0.0060,
    "HAARTadh": 0.36,
    "HAARTdisc": 0.023,
}
RaceClass2["ALL"] = {"Proportion": 0.0, "HAARTdisc": 0.018}

DemographicParams = {"WHITE": RaceClass1, "BLACK": RaceClass2}


"""
Partnership durations and
"""
sexualDurations = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
sexualDurations[1] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
sexualDurations[2] = {"p_value": (0.323 + 0.262 + 0.116), "min": 7, "max": 12}
sexualDurations[3] = {"p_value": (0.323 + 0.262 + 0.116 + 0.121), "min": 13, "max": 24}
sexualDurations[4] = {"p_value": (0.323 + 0.262 + 0.116 + 0.121 + 0.06), "min": 25, "max": 36}
sexualDurations[5] = {"min": 37, "max": 48}

needleDurations = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
needleDurations[1] = {"p_value": 1.0, "min": 1, "max": 6}
needleDurations[2] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[3] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[4] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[5] = {"min": 1, "max": 6}

PartnershipDurations = {"SEX": sexualDurations, "NEEDLE": needleDurations}


"""
Sexual and injection transmission probabilities
"""
SexTrans = {"MSM": {}, "HM": {}, "HF": {}}
SexTrans["MSM"] = {
    "0": 0.00745,
    "1": 0.00745,
    "2": 0.00745,
    "3": 0.00745,
    "4": 0.00745,
    "5": 0.000298,
}
SexTrans["HM"] = {"0": 0.001, "1": 0.001, "2": 0.0008, "3": 0.0004, "4": 0.0002, "5": 0.0001}
SexTrans["HF"] = {"0": 0.001, "1": 0.001, "2": 0.0008, "3": 0.0004, "4": 0.0002, "5": 0.0001}

NeedleTrans = {0: 0.007, "1": 0.007, "2": 0.0056, "3": 0.0028, "4": 0.0014, "5": 0.0002}

TransmissionProbabilities = {"SEX": SexTrans, "NEEDLE": NeedleTrans}
