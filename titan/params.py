__author__ = "MaximilianKing"

from typing import Sequence, List, Dict, Optional, Any
from copy import deepcopy

"""
Main model parameters.
"""

####################
PROCESSES = 1  # number of processes in parallel (quadcore)
rSeed = (
    0  # seed for random number generator (0 for pure random, -1 for stepwise up to N_NC
)
rSeed_pop = 0
rSeed_net = 0
rSeed_run = 0
N_MC = 1  # total number of iterations (Monte Carlo runs)
N_REPS = 1
N_POP = 100  # population size
TIME_RANGE = 12  # total time steps to iterate
burnDuration = 6  # total time for burning in period (equillibration)
model = "Custom"  # Model Type for fast flag toggling
setting = "Phil2005"
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
    "deathReport",
    "incarReport",
    "newlyhighriskReport",
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
cal_SexualActScaling = 0.45  # Sexual acts  scaling factor
cal_pXmissionScaling = 1.0  # 0.92 # Global transmission probability scaling factor
cal_AcuteScaling = 4.3  # Infectivity multiplier ratio for Acute status infections
cal_RR_Dx = 0.0  # Risk reduction in transmission probability for agents diagnosed
cal_RR_HAART = 1.0  # Scaling factor for effectiveness of ART therapy on xmission P
cal_TestFreq = 0.3  # Scaling factor for testing frequency
cal_Mortality = 0.5  # Scaling factor for all cause mortality rates
cal_ProgAIDS = 1.0  # Scaling factor for all progression to AIDS from HIV rates
cal_ART_cov = 0.70  # Scaling factor for enrollment on ART probability
cal_IncarP = 1.0  # Scaling factor for probability of becoming incarcerated
cal_raceXmission = (
    1.0  # Scaling factor for increased STI transmission P comparing race1/race2
)
cal_ptnrSampleDepth = 100  # Sampling depth for partnering algorithm.

"""
High risk params
"""
HR_partnerScale = 300  # Linear increase to partner number during HR period
HR_proportion = 0.3  # Proportion of people who enter HR group when partner incarcerated
HR_M_dur = 6  # Duration of high risk for males
HR_F_dur = 6  # Duration of high risk for females
condomUseType = "Race"  # Race or Acts
HIV_MSMW = 0.0

# Misc. params
flag_AssortativeMix = True
AssortMixType = "Race"
flag_AgeAssortMix = False
flag_RaceAssortMix = True
AssortMixCoeff = 0.75  # Proportion of race1 mixing with race2 when partnering.
safeNeedleExchangePrev = 1.0  # Prevalence scalar on SNE
initTreatment = 0
treatmentCov = 0.0
maxComponentSize = 1000
minComponentSize = 1

"""
Vaccine params
"""
vaccine_type = "RV144"
booster = True
vaccine_start = 1

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
inc_PtnrDissolution = 0.55
inc_treatment_startdate = 48  # Timestep where inc treatment can begin
inc_treatment_dur = (
    12  # Duration for which agents are forced on respective treatment post release
)
inc_treat_set = ["HM"]  # Set of agent classifiers effected by HR treatment
inc_treat_HRsex_beh = True  # Remove sexual higrisk behaviour during treatment duration
inc_treat_IDU_beh = True  # Remove IDU behav:iour during treatment duration
inc_treat_RIC = False  # Force retention in care of ART therapy

# PrEP params
PrEP_type = ["Oral", "Vaccine"]  # Oral/Inj PrEP and/or vaccine
PrEP_Target = (
    0.5  # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
)
PrEP_startT = 0  # Start date for PrEP program (0 for start of model)
PrEP_Adherence = 0.82  # Probability of being adherent
PrEP_AdhEffic = 0.96  # Efficacy of adherence PrEP
PrEP_NonAdhEffic = 0.76  # Efficacy of non-adherence PrEP
PrEP_falloutT = 0  # During PrEP remains effective post discontinuation
PrEP_resist = 0.01
PrEP_disc = 0.15
PrEP_target_model = (
    "CDCwomen"  # Allcomers, Clinical, Allcomers, HighPN5, HighPN10, SRIns, SR,Rec, MSM
)
PrEP_clinic_cat = "Mid"  # If clinical target model, which category does it follow

if "Oral" in PrEP_type:
    PrEP_Adherence = 1.0
    PrEP_AdhEffic = 0.96
    PrEP_NonAdhEffic = 0.76
    PrEP_falloutT = 1
    PrEP_disc = 0.15
    PrEP_peakLoad = 1.0
    PrEP_halflife = 1.0
elif "Inj" in PrEP_type:
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
    flag_high_risk = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_booster = False
elif model == "Incar":
    flag_incar = True
    flag_PrEP = False
    flag_high_risk = True
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_booster = False
elif model == "NoIncar":
    flag_incar = False
    flag_PrEP = False
    flag_high_risk = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_booster = False
elif model == "VaccinePrEP":
    flag_incar = False
    flag_PrEP = True
    flag_high_risk = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_booster = True
elif model == "Custom":
    flag_incar = False
    flag_PrEP = True
    flag_high_risk = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_booster = False
    flag_PCA = False
    init_with_vaccine = False

agentSexTypes = ["HM", "HF", "MSM", "WSW", "MTF"]
agentPopulations = deepcopy(agentSexTypes)
agentPopulations.append("IDU")

"""
RaceClass is a distinct racial/ethnic/social classification for demographics of the population.
ID is the specific mode of partnership the agent engages in (ie MSM, HM, HF, PWID)
RaceClass agent classifier template
"""
RC_template: Dict[str, Any] = {
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
    "HAARTprev": 0.0,
    "HAARTadh": 0.0,  # Adherence to ART therapy
    "HAARTdisc": 0.0,  # Probability of discontinuing ART therapy
    "EligSE_PartnerType": "MSM",  # List of agent SO types the agent cant partner with
    "PrEPdisc": 0.0,  # Probability of discontinuing PrEP treatment
    "HighRiskPrev": 0.0,
    "PrEPadh": 1.0,
    "PrEP_coverage": 0,
    "boosterInterval": 0,
    "boosterProb": 0,
    "vaccinePrev": 0,
    "vaccineInit": 0,
    "HighRiskDuration": 10,
}

RaceClass1: Dict[str, Any] = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}, "ALL": {}}
RaceClass2: Dict[str, Any] = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}, "ALL": {}}
for a in agentPopulations:
    RaceClass1[a] = dict(RC_template)
    RaceClass2[a] = dict(RC_template)

RaceClass2["HM"].update({"POP": 0.2})
RaceClass2["HF"].update({"POP": 0.8})

RaceClass1["HM"].update(
    {
        "POP": 0.0,
        "HIV": 0.0369,
        "AIDS": 0.6780,
        "HAARTprev": 0.41,
        "INCARprev": 0.0274,
        "TestedPrev": 0.90,
        "NUMPartn": 1.5,
        "NUMSexActs": 5.0,
        "SAFESEX": 0.89,
        "NEEDLESH": 0.43,
        "HIVTEST": 0.034,
        "INCAR": 0.001,
        "HAARTadh": 0.405,
        "HAARTdisc": 0.000,
        "PrEPdisc": 0.0000,
        "EligSE_PartnerType": "HF",
    }
)

RaceClass1["MSM"].update(
    {
        "POP": 1.00,
        "HIV": 0.5,
        "AIDS": 0.07,
        "HAARTprev": 0.583,
        "INCARprev": 0.000,
        "TestedPrev": 0.826,
        "mNPart": 7.0,
        "NUMPartn": 7.0,
        "NUMSexActs": 5.0,
        "SAFESEX": 0.432,
        "NEEDLESH": 0.43,
        "HIVTEST": 0.055,
        "INCAR": 0.00,  # 0.00014,
        "HAARTadh": 0.885,  # 0.693,#0.57,
        "HAARTdisc": 0.10,
        "PrEPdisc": 0.13,
        "EligSE_PartnerType": "MSM",
        "PrEPadh": 0.911,
        "PrEP_coverage": 0.0,
        "vaccinePrev": 1,
        "boosterInterval": 3,
        "boosterProb": 1.0,
    }
)

RaceClass1["ALL"].update(
    {"Proportion": 0.611, "HAARTdisc": 0.018, "PrEPdisc": 0.0, "AssortMixCoeff": 0.722}
)

# RaceClass2 = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}, 'ALL':{}}
RaceClass2["MSM"].update(
    {
        "POP": 0.00,  # 0.028,
        "HIV": 0.5,  # 0.434,
        "AIDS": 0.232,
        "HAARTprev": 0.627,
        "INCARprev": 0.00,
        "TestedPrev": 0.655,
        "mNPart": 5.0,
        "NUMPartn": 5.0,
        "NUMSexActs": 5.0,
        "SAFESEX": 0.312,
        "NEEDLESH": 0.27,
        "HIVTEST": 0.06,
        "INCAR": 0.00,  # 0.0011,
        "HAARTadh": 0.817,  # 0.598,#0.34,
        "HAARTdisc": 0.07,
        "PrEPdisc": 0.10,
        "EligSE_PartnerType": "MSM",
        "PrEPadh": 0.568,
        "PrEP_coverage": 0.0,
        "vaccinePrev": 1.0,
        "boosterInterval": 3,
        "boosterProb": 1.0,
    }
)
RaceClass2["HM"].update({"EligSE_PartnerType": "HF"})

RaceClass2["ALL"].update(
    {"Proportion": 0.389, "HAARTdisc": 0.018, "PrEPdisc": 0.0, "AssortMixCoeff": 0.765}
)

DemographicParams = {"WHITE": RaceClass1, "BLACK": RaceClass2}

"""
Partnership durations and
"""
sexualDurations: Dict[int, Any] = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
sexualDurations[1] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
sexualDurations[2] = {"p_value": (0.323 + 0.262 + 0.116), "min": 7, "max": 12}
sexualDurations[3] = {"p_value": (0.323 + 0.262 + 0.116 + 0.121), "min": 13, "max": 24}
sexualDurations[4] = {
    "p_value": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
sexualDurations[5] = {"min": 37, "max": 48}

needleDurations: Dict[int, Any] = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
needleDurations[1] = {"p_value": 1.0, "min": 1, "max": 6}
needleDurations[2] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[3] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[4] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[5] = {"min": 1, "max": 6}

PartnershipDurations = {"SEX": sexualDurations, "NEEDLE": needleDurations}

"""
Partnership acts and
"""
sexualFrequency: Dict[int, Any] = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
sexualFrequency[1] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
sexualFrequency[2] = {"p_value": (0.323 + 0.262 + 0.116), "min": 7, "max": 12}
sexualFrequency[3] = {"p_value": (0.323 + 0.262 + 0.116 + 0.121), "min": 13, "max": 24}
sexualFrequency[4] = {
    "p_value": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
sexualFrequency[5] = {
    "p_value": 1.0,
    "min": 37,
    "max": 48,
}

needleFrequency: Dict[int, Any] = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
needleFrequency[1] = {"p_value": 1.0, "min": 1, "max": 6}
needleFrequency[2] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleFrequency[3] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleFrequency[4] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleFrequency[5] = {"min": 1, "max": 6}

PartnershipDurations = {"SEX": sexualDurations, "NEEDLE": needleDurations}

"""
Sexual and injection transmission probabilities
"""
SexTrans: Dict[str, Any] = {"MSM": {}, "HM": {}, "HF": {}}

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

TransmissionProbabilities: Dict[str, Any] = {"SEX": SexTrans, "NEEDLE": NeedleTrans}

# ageMatrix[race][2]['Prop']

"""
Age bin distributions and HIV if utilized
"""
ageMatrix: Dict[str, Any] = {"WHITE": {}, "BLACK": {}}
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


"""
Age mixing matrix for assortative mixing by age
"""
mixingMatrix: Dict[int, Any] = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
mixingMatrix[1] = {1: 0.500, 2: 0.226, 3: 0.123, 4: 0.088, 5: 0.064}
mixingMatrix[2] = {1: 0.156, 2: 0.500, 3: 0.185, 4: 0.099, 5: 0.061}
mixingMatrix[3] = {1: 0.074, 2: 0.162, 3: 0.500, 4: 0.184, 5: 0.079}
mixingMatrix[4] = {1: 0.057, 2: 0.093, 3: 0.199, 4: 0.500, 5: 0.150}
mixingMatrix[5] = {1: 0.062, 2: 0.086, 3: 0.128, 4: 0.224, 5: 0.500}

"""
Clinic bins for targetting strategies
Bins represent partner numbers of the following category 0:0, 1:1, 2:2,  3:3-4, 4:5-9, 5:10+
"""
clinicAgents: Dict[str, Any] = {"Low": {}, "Mid": {}, "High": {}}
clinicAgents["Mid"] = {
    0: {"Prob": 0.0, "min": 0, "max": 0},
    1: {"Prob": 0.054, "min": 0, "max": 1},
    2: {"Prob": 0.061, "min": 1, "max": 2},
    3: {"Prob": 0.168, "min": 3, "max": 4},
    4: {"Prob": 0.246, "min": 5, "max": 9},
    5: {"Prob": 0.471, "min": 10, "max": 120},
}


"""
Bond Params
"""
nonSex = 0.308
multiplex = 0.105
sexualOnly = 0.587
bond_type = ["social"]
mean_partner_type = "bins"

"""
Peer change params
"""

LAI_chance = 0
attitude = {0: 0.167, 1: 0.082, 2: 0.184, 3: 0.139, 4: 0.429}
PCA_PrEP = 0.30 * (
    1 - 0.575
)  # chance of attempting * chance of initiating oral or inj PrEP
opinion_threshold = (
    3.0  # opinion needed to initiate PrEP on a 0-4 scale (translated from 1-5 scale)
)
pcaChoice = "bridge"  # eigenvector or bridge, how the PCA is selected
awarenessProb = 0.055  # static probability of becoming spontaneously aware of PrEP
starting_awareness = 0.00  # awareness of PrEP at t0
knowledgeTransmission = (
    0.01  # per-act probability of knowledge change in unaware partner
)
opinionTransmission = (
    0.005  # per-act probability of opinion change in less-prominent partner
)
interactionProb: Dict[str, Any] = {
    "sexOnly": {},
    "multiplex": {},
    "social": {},
}  # prob of interaction per timestep (or at relationship formation for sexual)
interactionProb["sexOnly"][1] = {"pvalue": 1.00, "min": 0, "max": 0}


interactionProb["multiplex"][1] = {"pvalue": 0.306, "min": 0, "max": 0}
interactionProb["multiplex"][2] = {"pvalue": 0.144, "min": 1, "max": 1}
interactionProb["multiplex"][3] = {"pvalue": 0.067, "min": 2, "max": 2}
interactionProb["multiplex"][4] = {"pvalue": 0.106, "min": 4, "max": 4}
interactionProb["multiplex"][5] = {"pvalue": 0.150, "min": 5, "max": 29}
interactionProb["multiplex"][6] = {"pvalue": 0.228, "min": 30, "max": 30}

interactionProb["social"][1] = {"pvalue": 0.253, "min": 0, "max": 0}
interactionProb["social"][2] = {"pvalue": 0.123, "min": 1, "max": 1}
interactionProb["social"][3] = {"pvalue": 0.060, "min": 2, "max": 2}
interactionProb["social"][4] = {"pvalue": 0.140, "min": 4, "max": 4}
interactionProb["social"][5] = {"pvalue": 0.168, "min": 5, "max": 29}
interactionProb["social"][6] = {"pvalue": 0.256, "min": 30, "max": 30}

partnerNumber: Dict[int, Any] = {
    0: 0.083,
    1: 0.181,
    2: 0.229,
    3: 0.172,
    4: 0.112,
    5: 0.102,
    6: 0.071,
    7: 0.028,
    8: 0.019,
    9: 0.005,
}
