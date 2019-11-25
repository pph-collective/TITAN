__author__ = "MaximilianKing"

"""
Main model parameters.
"""

####################
PROCESSES = 1  # number of processes in parallel (quadcore)
rSeed_pop = (
    0  # seed for random number generator (0 for pure random, -1 for stepwise up to N_NC
)
rSeed_net = 0
rSeed_run = 0
N_MC = 1  # total number of iterations (Monte Carlo runs)
N_REPS = 1
N_POP = 1000  # population size
TIME_RANGE = 10  # total time steps to iterate
burnDuration = 1
model = "PrEP"  # Model Type for fast flag toggling
setting = "AtlantaMSM"
network_type = "scale_free"
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
drawFigures = False
calcComponentStats = False
flag_agentZero = False

reports = [
    "incidenceReport",
    "prevalenceReport",
    "deathReport",
    "incarReport",
    "iduReport",
    "highriskReport",
    "newlyhighriskReport",
    "sexReport",
    "prepReport",
    "basicReport",
]

"""
Calibration scaling parameters for fitting to empirical data
"""

PARTNERTURNOVER = (
    1.0 / 7.5
)  # Partner acquisition parameters (higher number more partnering)
cal_NeedlePartScaling = 1.00  # IDU partner number scaling
cal_NeedleActScaling = 1.00  # IDU act frequency scaling factor
cal_SexualPartScaling = 1.0  # Sexual partner number scaling factor
cal_SexualActScaling = 2.0  # Sexual acts  scaling factor
cal_pXmissionScaling = 1.0  # 0.92 # Global transmission probability scaling factor
cal_AcuteScaling = 4.3  # Infectivity multiplier ratio for Acute status infections
cal_RR_Dx = 0.0  # Risk reduction in transmission probability for agents diagnosed
cal_RR_HAART = 1.0  # Scaling factor for effectiveness of ART therapy on xmission P
cal_TestFreq = 0.3  # Scaling factor for testing frequency
cal_Mortality = 0.5  # Scaling factor for all cause mortality rates
cal_ProgAIDS = 0.05  # Scaling factor for all progression to AIDS from HIV rates
cal_ART_cov = 0.4  # Scaling factor for enrollment on ART probability
cal_IncarP = 1.0
cal_raceXmission = 4.0
cal_ptnrSampleDepth = 100

# High risk params
HR_partnerScale = 300  # Linear increase to partner number during HR period
HR_proportion = 0.3  # Proportion of people who enter HR group when partner incarcerated
HR_M_dur = 6  # Duration of high risk for males
HR_F_dur = 6  # Duration of high risk for females
condomUseType = "Race"  # Race or Acts

# Misc. params
flag_AssortativeMix = True
AssortMixType = "Race"
flag_AgeAssortMix = False
flag_RaceAssortMix = True
AssortMixCoeff = 0.75  # Proportion of race1 mixing with race2 when partnering.
safeNeedleExchangePrev = 1.0  # Prevalence scalar on SNE
initTreatment = 999999
treatmentCov = 0.0

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
PrEP_type = "Oral"  # Oral/Inj PrEP modes
PrEP_Target = (
    0.00  # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
)
PrEP_startT = 0  # Start date for PrEP program (0 for start of model)
PrEP_Adherence = 0.82  # Probability of being adherent
PrEP_AdhEffic = 0.96  # Efficacy of adherence PrEP
PrEP_NonAdhEffic = 0.76  # Efficacy of non-adherence PrEP
PrEP_falloutT = 0  # During PrEP remains effective post discontinuation
PrEP_resist = 0.01
PrEP_disc = 0.15
PrEP_target_model = (
    "Allcomers"  # Clinical, Allcomers, HighPN5, HighPN10, SRIns, SR,CDC,Racial
)
PrEP_init_var1 = 0.05
PrEP_init_var2 = 0.025
PrEP_clinic_cat = "Racial"

if PrEP_type == "Oral":
    PrEP_Adherence = "byRace"
    PrEP_AdhEffic = 0.96
    PrEP_NonAdhEffic = 0.76
    PrEP_falloutT = 1
    PrEP_disc = 0.15
elif PrEP_type == "Inj":
    PrEP_Adherence = 1.0
    PrEP_AdhEffic = 1.0
    PrEP_NonAdhEffic = 1.00
    PrEP_falloutT = 12
    PrEP_disc = 0.04
    PrEP_peakLoad = 4.91
    PrEP_halflife = 40.0

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
    flag_staticN = False
elif model == "Incar":
    flag_incar = True
    flag_PrEP = False
    flag_HR = True
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
elif model == "NoIncar":
    flag_incar = False
    flag_PrEP = False
    flag_HR = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
elif model == "Custom":
    flag_incar = False
    flag_PrEP = True
    flag_HR = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False

agentPopulations = ["MSM"]
agentSexTypes = ["MSM"]

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
    "HAARTprev": 0.0,
    "HAARTadh": 0.0,  # Adherence to ART therapy
    "HAARTdisc": 0.0,  # Probability of discontinuing ART therapy
    "EligSE_PartnerType": None,  # List of agent SO types the agent cant partner with
    "PrEPdisc": 0.0,  # Probability of discontinuing PrEP treatment
    "HighRiskPrev": 0.0,
    "PrEPadh": 1.0,
}

# RaceClass1 = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}, 'ALL':{}}
RaceClass1 = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}, "ALL": {}}
RaceClass2 = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}, "ALL": {}}
for a in ["MSM", "HM", "HF", "IDU"]:
    RaceClass1[a] = dict(RC_template)
    RaceClass2[a] = dict(RC_template)

RaceClass1["MSM"]["POP"] = 1.0
RaceClass1["MSM"]["HIV"] = 0.4
# StratW['MSM'] = {'POP':0.035, 'HIV':0.132, 'AIDS':0.048, 'HAARTprev':0.57, 'INCARprev':0.005, 'TestedPrev':0.84}

RaceClass1["MSM"].update(
    {
        "POP": 1.00,
        "HIV": 0.132,
        "AIDS": 0.07,
        "HAARTprev": 0.410,  # 0.895,
        "INCARprev": 0.000,
        "TestedPrev": 0.826,
        "mNPart": 7.0,
        "NUMPartn": 7.0,
        "NUMSexActs": 5.0,
        "UNSAFESEX": 0.432,
        "NEEDLESH": 0.43,
        "HIVTEST": 0.055,
        "INCAR": 0.00,  # 0.00014,
        "HAARTadh": 0.885,  # 0.693,#0.57,
        "HAARTdisc": 0.008,
        "PrEPdisc": 0.13,
        "EligSE_PartnerType": "MSM",
        "PrEPadh": 0.911,
    }
)

RaceClass1["ALL"].update(
    {"Proportion": 0.611, "HAARTdisc": 0.018, "PrEPdisc": 0.0, "AssortMixCoeff": 0.722}
)

# RaceClass2 = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}, 'ALL':{}}
RaceClass2["MSM"].update(
    {
        "POP": 1.00,  # 0.028,
        "HIV": 0.434,
        "AIDS": 0.232,
        "HAARTprev": 0.309,  # 0.845,
        "INCARprev": 0.00,
        "TestedPrev": 0.655,
        "mNPart": 5.0,
        "NUMPartn": 5.0,
        "NUMSexActs": 5.0,
        "UNSAFESEX": 0.312,
        "NEEDLESH": 0.27,
        "HIVTEST": 0.06,
        "INCAR": 0.00,  # 0.0011,
        "HAARTadh": 0.817,  # 0.598,#0.34,
        "HAARTdisc": 0.01,
        "PrEPdisc": 0.15,
        "EligSE_PartnerType": "MSM",
        "PrEPadh": 0.568,
    }
)

RaceClass2["ALL"].update(
    {"Proportion": 0.389, "HAARTdisc": 0.018, "PrEPdisc": 0.0, "AssortMixCoeff": 0.765}
)

DemographicParams = {"WHITE": RaceClass1, "BLACK": RaceClass2}

"""
Partnership durations and
"""
sexualDurations = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
sexualDurations[1] = {"p_value": 0.456, "min": 1, "max": 3}
sexualDurations[2] = {"p_value": (0.456 + 0.300), "min": 3, "max": 12}
sexualDurations[3] = {"p_value": (0.456 + 0.300 + 0.245), "min": 13, "max": 24}
sexualDurations[4] = {"p_value": (0.281 + 0.209 + 0.281 + 0.230), "min": 13, "max": 24}
sexualDurations[5] = {"min": 13, "max": 24}

needleDurations = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
needleDurations[1] = {"p_value": 1.0, "min": 1, "max": 6}
needleDurations[2] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[3] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[4] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[5] = {"min": 1, "max": 6}

PartnershipDurations = {"SEX": sexualDurations, "NEEDLE": needleDurations}

"""
Partnership acts and
"""
# todo FINISH THESE AND IMPORT HTEM INTO SEXACTS
sexualFrequency = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
sexualFrequency[1] = {"p_value": (0.244), "min": 1, "max": 2}
sexualFrequency[2] = {"p_value": (0.244 + 0.493), "min": 3, "max": 4}
sexualFrequency[3] = {"p_value": (0.244 + 0.493 + 0.123), "min": 5, "max": 12}
sexualFrequency[4] = {"p_value": (0.244 + 0.493 + 0.123 + 0.14), "min": 13, "max": 20}
sexualFrequency[5] = {"p_value": (0.244 + 0.493 + 0.123 + 0.14), "min": 13, "max": 20}
sexualFrequency[6] = {"p_value": (0.244 + 0.493 + 0.123 + 0.14), "min": 13, "max": 20}
sexualFrequency[7] = {"p_value": (0.244 + 0.493 + 0.123 + 0.14), "min": 13, "max": 20}
sexualFrequency[8] = {"p_value": (0.244 + 0.493 + 0.123 + 0.14), "min": 13, "max": 20}
sexualFrequency[9] = {"min": 13, "max": 20}

needleFrequency = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
needleFrequency[1] = {"p_value": 1.0, "min": 1, "max": 6}
needleFrequency[2] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleFrequency[3] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleFrequency[4] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleFrequency[5] = {"min": 1, "max": 6}

PartnershipDurations = {"SEX": sexualDurations, "NEEDLE": needleDurations}

"""
Sexual and injection transmission probabilities
"""
SexTrans = {"MSM": {}, "HM": {}, "HF": {}}
SexTrans["MSM"] = {
    "0": 0.00745,
    "1": 0.00745 * 0.81,
    "2": 0.00745 * 0.81,
    "3": 0.00745 * 0.81,
    "4": 0.00745 * 0.81,
    "5": 0.000,
}
SexTrans["HM"] = {
    "0": 0.001,
    "1": 0.001,
    "2": 0.0008,
    "3": 0.0004,
    "4": 0.0002,
    "5": 0.000,
}
SexTrans["HF"] = {
    "0": 0.001,
    "1": 0.001,
    "2": 0.0008,
    "3": 0.0004,
    "4": 0.0002,
    "5": 0.000,
}

NeedleTrans = {
    "0": 0.007,
    "1": 0.007,
    "2": 0.0056,
    "3": 0.0028,
    "4": 0.0014,
    "5": 0.0002,
}

TransmissionProbabilities = {"SEX": SexTrans, "NEEDLE": NeedleTrans}

# ageMatrix[race][2]['Prop']

ageMatrix = {"WHITE": {}, "BLACK": {}}
ageMatrix["WHITE"] = {
    "Prop": {
        0: 0.0,
        1: 0.085,
        2: 0.085 + 0.206,
        3: 0.085 + 0.206 + 0.222,
        4: 0.085 + 0.206 + 0.222 + 0.493,
        5: 2,
    },
    "HIV": {0: 0.0, 1: 0.006, 2: 0.029, 3: 0.055, 4: 0.069, 5: 0.025},
}
ageMatrix["BLACK"] = {
    "Prop": {
        0: 0.0,
        1: 0.105,
        2: 0.105 + 0.225,
        3: 0.105 + 0.225 + 0.215,
        4: 0.105 + 0.225 + 0.215 + 0.455,
        5: 0.28 + 0.24 + 0.19 + 0.15 + 0.14,
    },
    "HIV": {0: 0.0, 1: 0.144, 2: 0.144, 3: 0.250, 4: 0.377, 5: 0.194},
}

mixingMatrix = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
mixingMatrix[1] = {1: 0.500, 2: 0.226, 3: 0.123, 4: 0.088, 5: 0.064}
mixingMatrix[2] = {1: 0.156, 2: 0.500, 3: 0.185, 4: 0.099, 5: 0.061}
mixingMatrix[3] = {1: 0.074, 2: 0.162, 3: 0.500, 4: 0.184, 5: 0.079}
mixingMatrix[4] = {1: 0.057, 2: 0.093, 3: 0.199, 4: 0.500, 5: 0.150}
mixingMatrix[5] = {1: 0.062, 2: 0.086, 3: 0.128, 4: 0.224, 5: 0.500}

# Bins represent partner numbers of the following category 0:0, 1:1, 2:2,  3:3-4, 4:5-9, 5:10+
clinicAgents = {"Racial": {}, "Mid": {}, "High": {}}
clinicAgents["Racial"] = {
    0: {"Prob": 10.0, "min": 0, "max": 0},
    1: {"Prob": 10.054, "min": 0, "max": 1},
    2: {"Prob": 10.061, "min": 1, "max": 2},
    3: {"Prob": 10.168, "min": 3, "max": 4},
    4: {"Prob": 10.246, "min": 5, "max": 9},
    5: {"Prob": 10.471, "min": 10, "max": 120},
}
clinicAgents["Mid"] = {
    0: {"Prob": 0.0, "min": 0, "max": 0},
    1: {"Prob": 0.054, "min": 0, "max": 1},
    2: {"Prob": 0.061, "min": 1, "max": 2},
    3: {"Prob": 0.168, "min": 3, "max": 4},
    4: {"Prob": 0.246, "min": 5, "max": 9},
    5: {"Prob": 0.471, "min": 10, "max": 120},
}
