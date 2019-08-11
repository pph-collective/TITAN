__author__ = "MaximilianKing"


"""
Main model parameters.
"""

####################
PROCESSES = 1  # number of processes in parallel (quadcore)
rSeed = 1  # seed for random number generator (0 for pure random, -1 for stepwise up to N_NC
N_MC = 1  # total number of iterations (Monte Carlo runs)
N_POP = 30000  # population size
TIME_RANGE = 120  # total time steps to iterate
model = "Custom"  # Model Type for fast flag toggling
####################

"""
Output flags and settings
"""
outputDir = ""

startAgentList = False
endingAgentList = False
intermAgentList = False
intermPrintFreq = 10
MSMreport = True
HMreport = False
HFreport = False


"""
Calibration scaling parameters for fitting to empirical data
"""

PARTNERTURNOVER = 7.5  # Partner acquisition parameters (higher number more partnering)

cal_NeedleScaling = 1.0  # IDU transmission probability scaling factor
cal_SexualScaling = 1.0  # Sexual transmission probability scaling factor
cal_pXmissionScaling = 1.0  # 0.92 # Global transmission probability scaling factor
cal_AcuteScaling = 4.3  # Infectivity multiplier ratio for Acute status infections
cal_RR_Dx = 0.25  # Risk reduction in transmission probability for agents diagnosed
cal_RR_HAART = 3.0  # Scaling factor for effectiveness of ART therapy on xmission P
cal_TestFreq = 0.5  # Scaling factor for testing frequency
cal_Mortality = 0.05  # Scaling factor for all cause mortality rates
cal_ProgAIDS = 0.05  # Scaling factor for all progression to AIDS from HIV rates
cal_ART_cov = 1.45  # Scaling factor for enrollment on ART probability
cal_IncarP = 1.0
cal_raceXmission = 2.9

# High risk params
HR_partnerScale = 300  # Linear increase to partner number during HR period
HR_proportion = 0.3  # Proportion of people who enter HR group when partner incarcerated
HR_M_dur = 6  # Duration of high risk for males
HR_F_dur = 6  # Duration of high risk for females


# Misc. params
AssortMixCoeff = 0.05  # Proportion of race1 mixing with race2 when partnering.
safeNeedleExchangePrev = 1.0

# Incarceration params
inc_JailMax = 9
inc_JailMin = 1
inc_JailTestProb = 0.69
inc_PrisMax = 60
inc_PrisMin = 6
inc_PrisTestProb = 0.69
inc_PropPrison = 0.5
inc_ARTenroll = 0.51
inc_ARTadh = 0.21
inc_ARTdisc = 0.12
inc_Recidivism = 0.267

# PrEP params
PrEP_type = "Inj"  # Oral/Inj PrEP modes
PrEP_Target = 0.0  # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
PrEP_startT = 0  # Start date for PrEP program (0 for start of model)
PrEP_Adherence = 0.82  # Probability of being adherent
PrEP_AdhEffic = 0.96  # Efficacy of adherence PrEP
PrEP_NonAdhEffic = 0.76  # Efficacy of non-adherence PrEP
PrEP_falloutT = 0  # During PrEP remains effective post discontinuation
PrEP_resist = 0.0  # 1
PrEP_disc = 0.15

if PrEP_type == "Oral":
    PrEP_Adherence = 0.92
    PrEP_AdhEffic = 0.96
    PrEP_NonAdhEffic = 0.76
    PrEP_falloutT = 1
    PrEP_disc = 0.15
elif PrEP_type == "Inj":
    PrEP_Adherence = 1.0
    PrEP_AdhEffic = 0.90
    PrEP_NonAdhEffic = 1.00
    PrEP_falloutT = 9
    PrEP_disc = 0.04

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
elif model == "Custom":
    flag_incar = False
    flag_PrEP = True
    flag_HR = False
    flag_ART = True
    flag_DandR = True


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
    "UNSAFESEX": 0.0,  # Probability of engaging in unsafe sex (per act)
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
    "POP": 1.00,  # 0.035,
    "HIV": 0.132,
    "AIDS": 0.048,
    "HAARTprev": 0.57,
    "INCARprev": 0.005,
    "TestedPrev": 0.84,
    "mNPart": 3,
    "NUMPartn": 1.5,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.43,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.055,
    "INCAR": 0.00,  # 0.00014,
    "HAARTadh": 1.0,  # 0.57,
    "HAARTdisc": 0.008,
    "PrEPdisc": PrEP_disc,
}

# StratW['HM'] = {'POP':0.465, 'HIV':0.003, 'AIDS':0.001, 'HAARTprev':0.56, 'INCARprev':0.005, 'TestedPrev':0.81}
RaceClass1["HM"] = {
    "POP": 0.00,  # 0.465,
    "HIV": 0.003,
    "AIDS": 0.001,
    "HAARTprev": 0.56,
    "INCARprev": 0.005,
    "TestedPrev": 0.81,
    "mNPart": 3,
    "NUMPartn": 1.2,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.83,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.013,
    "INCAR": 0.00014,
    "HAARTadh": 0.57,
    "HAARTdisc": 0.01,
    "PrEPdisc": 0.0,
}

# StratW['HF'] = {'POP':0.50, 'HIV':0.002, 'AIDS':0.001, 'HAARTprev':0.53, 'INCARprev':0.0002, 'TestedPrev':0.86}
RaceClass1["HF"] = {
    "POP": 0.00,  # 0.50,
    "HIV": 0.002,
    "AIDS": 0.001,
    "HAARTprev": 0.53,
    "INCARprev": 0.0002,
    "TestedPrev": 0.86,
    "mNPart": 3,
    "NUMPartn": 1.2,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.83,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.015,
    "INCAR": 0.00003,
    "HAARTadh": 0.53,
    "HAARTdisc": 0.013,
    "PrEPdisc": 0.0,
}

# StratW['PWID'] = {'POP':0.005, 'HIV':0.19, 'AIDS':0.075, 'HAARTprev':0.52, 'INCARprev':0.019, 'TestedPrev':0.93}
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
    "UNSAFESEX": 0.72,
    "NEEDLESH": 0.43,
    "HIVTEST": 0.037,
    "INCAR": 0.00,  # 0.0050,
    "HAARTadh": 0.43,
    "HAARTdisc": 0.018,
    "PrEPdisc": 0.0,
}
RaceClass1["ALL"] = {"Proportion": 0.571, "HAARTdisc": 0.018, "PrEPdisc": 0.0}

RaceClass2 = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
# StratB['MSM'] = {'POP':0.028, 'HIV':0.432 , 'AIDS':0.048, 'HAARTprev':0.47, 'INCARprev':0.018, 'TestedPrev':0.84}
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
    "UNSAFESEX": 0.43,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.06,
    "INCAR": 0.00,  # 0.0011,
    "HAARTadh": 1.0,  # 0.34,
    "HAARTdisc": 0.01,
    "PrEPdisc": PrEP_disc,
}
# StratB['HM'] = {'POP':0.472, 'HIV':0.012 , 'AIDS':0.003, 'HAARTprev':0.44, 'INCARprev':0.018, 'TestedPrev':0.84}
RaceClass2["HM"] = {
    "POP": 0.00,  # 0.472,
    "HIV": 0.012,
    "AIDS": 0.003,
    "HAARTprev": 0.44,
    "INCARprev": 0.018,
    "TestedPrev": 0.84,
    "mNPart": 3,
    "NUMPartn": 3.5,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.70,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.023,
    "INCAR": 0.0011,
    "HAARTadh": 0.34,
    "HAARTdisc": 0.013,
    "PrEPdisc": 0.0,
}
# StratB['HF'] = {'POP':0.50, 'HIV':0.021 , 'AIDS':0.0009, 'HAARTprev':0.52, 'INCARprev':0.0006, 'TestedPrev':0.86}
RaceClass2["HF"] = {
    "POP": 0.00,  # 0.50,
    "HIV": 0.021,
    "AIDS": 0.0009,
    "HAARTprev": 0.52,
    "INCARprev": 0.0006,
    "TestedPrev": 0.86,
    "mNPart": 3,
    "NUMPartn": 3.5,
    "NUMSexActs": 5.0,
    "UNSAFESEX": 0.70,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.033,
    "INCAR": 0.00016,
    "HAARTadh": 0.50,
    "HAARTdisc": 0.017,
    "PrEPdisc": 0.0,
}

# StratB['PWID'] = {'POP':0.005, 'HIV':0.19 , 'AIDS':0.075, 'HAARTprev':0.44, 'INCARprev':0.063, 'TestedPrev':0.93}
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
    "UNSAFESEX": 0.69,
    "NEEDLESH": 0.27,
    "HIVTEST": 0.043,
    "INCAR": 0.00,  # 0.0060,
    "HAARTadh": 0.36,
    "HAARTdisc": 0.023,
    "PrEPdisc": 0.0,
}
RaceClass2["ALL"] = {"Proportion": 0.429, "HAARTdisc": 0.018, "PrEPdisc": 0.0}

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
SexTrans["MSM"] = {"0": 0.005, "1": 0.005, "2": 0.004, "3": 0.002, "4": 0.001, "5": 0.0001}
SexTrans["HM"] = {"0": 0.001, "1": 0.001, "2": 0.0008, "3": 0.0004, "4": 0.0002, "5": 0.0001}
SexTrans["HF"] = {"0": 0.001, "1": 0.001, "2": 0.0008, "3": 0.0004, "4": 0.0002, "5": 0.0001}

NeedleTrans = {"0": 0.007, "1": 0.007, "2": 0.0056, "3": 0.0028, "4": 0.0014, "5": 0.0002}

TransmissionProbabilities = {"SEX": SexTrans, "NEEDLE": NeedleTrans}
