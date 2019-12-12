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
TIME_RANGE = 120  # total time steps to iterate
burnDuration = 36
model = "StaticZero"  # Model Type for fast flag toggling
network_type = "scale_free"  # scale_free or max_k_comp_size
setting = "Scott"
####################


"""
Output flags and settings
"""
outputDir = ""

printIncidenceEvents = False  # this needs a description
printStartAgentList = False  # prints the full agent list at the beginning of model run, after burn TODO check this
printEndingAgentList = False  # prints the full agent list at the end of model run
printIntermAgentList = False  # prints the full agent list at intervals during model run, determined by interimPrintFreq
intermPrintFreq = (
    60  # if printing agent lists during model run, determines frequency of printouts
)
calcNetworkStats = False  # TODO: ask Aly abt this?
calcComponentStats = (
    False  # prints individual component starting HIV, PrEP coverage, ending HIV, etc.
)
# TODO: only ashley?
drawFigures = False  # TODO: check. Used?
drawEdgeList = True  # prints the full network edge list (TODO: at time?)
drawFigureColor = "MSW"  # TODO: what?

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

PARTNERTURNOVER = 0.2  # Partner acquisition parameters (higher number more partnering)
cal_NeedlePartScaling = 1.0  # Scaling IDU partner number
cal_NeedleActScaling = 2.0  # Scaling factor IDU act frequency
cal_SexualPartScaling = 1.0  # Scaling factor for sexual partner number
cal_SexualActScaling = 1.0  # Scaling factor for sexual acts
cal_pXmissionScaling = (
    1.0  # Global transmission probability scaling factor. Do not change for calibration
)
cal_AcuteScaling = 4.3  # Infectivity multiplier ratio for Acute status infections. REVIEW: WILL, is this set in stone?
cal_RR_Dx = 0.53  # Scaling factor for risk reduction in transmission probability for agents diagnosed
cal_RR_HAART = 1.0  # Scaling factor for effectiveness of ART therapy on transmission
cal_TestFreq = 1.0  # Scaling factor for testing/diagnosis frequency
cal_Mortality = 0.5  # Scaling factor for all cause mortality rates
cal_ProgAIDS = 1.0  # Scaling factor for all rates of progression to AIDS from HIV
cal_ART_cov = 1.0  # Scaling factor for probability of enrollment on ART
cal_IncarP = 1.0  # Scaling factor for incarceration probability
cal_raceXmission = 1.0  # Scaling factor for racial disparities in transmission not accounted for in other factors
cal_ptnrSampleDepth = 100

"""
High risk params
"""
HR_partnerScale = 300  # Linear increase to partner number during HR period
HR_proportion = 0.3  # Proportion of people who enter HR group when partner incarcerated
HR_M_dur = 6  # Duration of high risk for males
HR_F_dur = 6  # Duration of high risk for females

"""
Misc. params
"""

flag_AssortativeMix = False  # Boolean for if assortative mixing occurs at all
AssortMixType = None  # Other assortative mixing types
flag_AgeAssortMix = False  # Assortative mix by age
flag_RaceAssortMix = False  # Assortative mix by race
AssortMixCoeff = 0.8  # Proportion of following given assort mix rules
safeNeedleExchangePrev = 1.0  # Prevalence scalar on SNE
initTreatment = 10  # Requirement to start treatment
treatmentCov = 0.60  # Prop that receive treatment
limitComponentSize = False  # Maintains a maximum component size
maxComponentSize = 100  # Component size maximum if limited
minComponentSize = 2  # Doesn't work?
condomUseType = "Acts"  # Racial or Acts. Acts is standard. Determines how likelihood of condom use is decided
"""
Vaccine params
"""
vaccine_type = "RV144"  # RV144 or "HVTN702". Type of vaccine administered.
booster = (
    True  # Allows boosting for the vaccine a defined time from first administration
)
vaccine_start = (
    1  # Time that vaccine is started. Agents cannot be given the first dose of
)
# vaccine in any other time step.

"""
Incarceration params
"""
# TODO: is this even how this is done anymore?
inc_JailMax = 22
inc_JailMin = 8
inc_JailTestProb = 0.69
inc_PrisMax = 96
inc_PrisMin = 45
inc_PrisTestProb = 0.69
inc_PropPrison = 0.5
inc_ARTenroll = 0.51
inc_ARTadh = 0.21
inc_ARTdisc = 0.12
inc_Recidivism = 0.267
inc_PtnrDissolution = 0.55
inc_treatment_dur = (
    6  # Duration for which agents are forced on respective treatment post release
)
inc_treat_set = ["HM"]  # Set of agent classifiers effected by HR treatment
inc_treat_behavior = True  # Remove IDU behaviour during treatment duration
inc_treat_RIC = False  # Force retention in care of ART therapy

"""
PrEP params
"""
PrEP_type = [
    "Oral"
]  # Type(s) of PrEP/prevention measure included. Current: Oral, Inj, Vaccine
PrEP_Target = (
    0.000  # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
)
PrEP_startT = 0  # Start date for PrEP program (0 for start of model)
PrEP_Adherence = 0.82  # Probability of being fully adherent to PrEP
PrEP_AdhEffic = 0.96  # Efficacy of PrEP with full adherence
PrEP_NonAdhEffic = 0.76  # Efficacy of PrEP with partial adherence
PrEP_falloutT = 0  # Duration PrEP remains effective after discontinuation
PrEP_resist = (
    0.01  # PrEP resistance probability  TODO: what exactly is this and how does it work
)
PrEP_disc = 0.15  # Probability of agent discontinuing PrEP in a given time step
PrEP_target_model = (
    "Allcomers"  # Clinical, Allcomers, HighPN5, HighPN10, SRIns, SR, Rec
)
PrEP_clinic_cat = "Mid"
# TODO the above get overwritten by the below?
if PrEP_type == "Oral":
    PrEP_Adherence = 0.923
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
if model == "PrEP":
    flag_incar = False
    flag_PrEP = True
    flag_HR = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_agentZero = False
    flag_booster = False

elif model == "Incar":
    flag_incar = True
    flag_PrEP = False
    flag_HR = True
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_agentZero = False
    flag_booster = False

elif model == "NoIncar":
    flag_incar = False
    flag_PrEP = False
    flag_HR = True
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_agentZero = False
    flag_booster = False

elif model == "StaticZero":
    flag_incar = False
    flag_PrEP = False
    flag_HR = False
    flag_ART = False
    flag_DandR = False
    flag_staticN = True
    flag_agentZero = True
    flag_booster = False

elif model == "Custom":
    flag_incar = False
    flag_PrEP = False
    flag_HR = False
    flag_ART = False
    flag_DandR = False
    flag_staticN = True
    flag_agentZero = False
    flag_booster = False

agentPopulations = [
    "HM",
    "HF",
]  # Populations in the model, possibilities HM, HF, MSM, WSW, MTF
agentSexTypes = ["HM", "HF", "MSM"]  # redundant w above?
"""
RaceClass is a distinct racial/ethnic/social classification for demographics of the population.
ID is the specific mode of partnership the agent engages in (ie MSM, HM, HF, PWID)
RaceClass agent classifier template
"""
RC_template = {
    "Race": None,  # Race of demographic
    "Class": None,  # Classification of networking
    "POP": 0.0,  # Percentage of total raceclass population that are ID
    "HIV": 0.0,  # Proportion of total ID population that are HIV
    "AIDS": 0.0,  # Proportion of total HIV_ID that are AIDS
    "HAARTprev": 0.0,  # Proportion of HIV_TESTED_ID that are enrolled on ART
    "INCARprev": 0.0,  # Proportion of ID that are incarcerated
    "TestedPrev": 0.0,  # Proportion of HIV_ID that are tested positively
    "HighRiskPrev": 0.0,  # Proportion of agents that are HR
    "mNPart": 0.0,  # Mean number of sex partners
    "NUMPartn": 0.0,  # Number of partners (redundant)
    "NUMSexActs": 0.0,  # Mean number of sex acts with each partner
    "UNSAFESEX": 0.0,  # Probability of engaging in unsafe sex (per act)
    "NEEDLESH": 0.0,  # Probability of sharing syringes during join drug use (per act)
    "HIVTEST": 0.0,  # Probability of testing for HIV
    "INCAR": 0.0,  # Probability of becoming incarcerated (rate)
    "HAARTadh": 0.0,  # Adherence to ART therapy
    "HAARTdisc": 0.0,  # Probability of discontinuing ART therapy
    "PrEPprev": 0.0,  # Proportion of HIV- that are enrolled on PrEP
    "PrEPdisc": 0.0,  # Probability of discontinuing PrEP treatment
    "EligSE_PartnerType": None,  # List of agent SO types the agent cant partner with
    "AssortMixMatrix": [],  # List of assortMix Matrix to be zipped with EligPart
    "boosterInterval": 9999999,  # Interval for booster in time since last vaccination
    "boosterProb": 0,  # Probability of eligible agent getting a vaccine booster
    "vaccinePrev": 0,  # Prevalence of vaccine in the population at vaccine start time
}

RC_allTemplate = {
    "Proportion": 1.00,  # Proportion of total population that is raceclass
    "HAARTdisc": 0.018,  # Overall HAART discontinuation probability
    "PrEPdisc": 0.0,  # Overall PrEP discontinuation probability
    "AssortMixCoeff": 1.0,  # Proportion RC mixes with other raceclass
}

RaceClass1 = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}, "ALL": {}}
RaceClass2 = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}, "ALL": {}}
for a in ["MSM", "HM", "HF", "IDU"]:
    RaceClass1[a] = dict(RC_template)
    RaceClass2[a] = dict(RC_template)

RaceClass1["HM"].update(
    {
        "POP": 0.49,
        "HIV": 0.0014,
        "AIDS": 0.6780,
        "HAARTprev": 0.41,
        "INCARprev": 0.0274,
        "TestedPrev": 0.90,
        "HighRiskPrev": 0.0,
        "NUMPartn": 1.5,
        "NUMSexActs": 13.4,
        "UNSAFESEX": 0.89,
        "NEEDLESH": 0.43,
        "HIVTEST": 0.034,
        "INCAR": 0.001,
        "HAARTadh": 0.405,
        "HAARTdisc": 0.000,
        "PrEPdisc": 0.0000,
        "EligSE_PartnerType": "HF",
    }
)

RaceClass1["HF"].update(
    {
        "POP": 0.51,
        "HIV": 0.0004,
        "AIDS": 0.573,
        "HAARTprev": 0.47,
        "INCARprev": 0.000,
        "TestedPrev": 0.90,
        "HighRiskPrev": 0.0,
        "NUMPartn": 0.5,
        "NUMSexActs": 12.74,
        "UNSAFESEX": 0.43,
        "NEEDLESH": 0.43,
        "HIVTEST": 0.034,
        "INCAR": 0.00,
        "HAARTadh": 0.405,
        "HAARTdisc": 0.000,
        "PrEPdisc": PrEP_disc,
        "EligSE_PartnerType": "HM",
    }
)

RaceClass1["IDU"].update(
    {
        "POP": 0.017,
        "HIV": 0.000,
        "AIDS": 0.6780,
        "HAARTprev": 0.41,
        "INCARprev": 0.0274,
        "TestedPrev": 0.90,
        "HighRiskPrev": 0.0,
        "NUMPartn": 0.5,
        "NUMSexActs": 5.0,
        "UNSAFESEX": 0.89,
        "NEEDLESH": 0.63,
        "HIVTEST": 0.055,
        "INCAR": 0.001,
        "HAARTadh": 0.405,
        "HAARTdisc": 0.000,
        "PrEPadh": 0.55,
        "PrEPdisc": 0.0000,
        "EligSE_PartnerType": "IDU",
    }
)

RaceClass1["ALL"].update(
    {"Proportion": 1.00, "HAARTdisc": 0.018, "PrEPdisc": 0.0, "AssortMixCoeff": 1.0}
)

RaceClass2["ALL"].update(
    {"Proportion": 0.00, "HAARTdisc": 0.018, "PrEPdisc": 0.0, "AssortMixCoeff": 1.0}
)

DemographicParams = {"WHITE": RaceClass1, "BLACK": RaceClass2}


"""
Partnership duration distribution bins
"""
# REVIEW: is this current logic?
sexualDurations = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
sexualDurations[1] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
sexualDurations[2] = {"p_value": (0.323 + 0.262 + 0.116), "min": 7, "max": 12}
sexualDurations[3] = {"p_value": (0.323 + 0.262 + 0.116 + 0.121), "min": 13, "max": 24}
sexualDurations[4] = {
    "p_value": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
sexualDurations[5] = {"min": 37, "max": 48}

needleDurations = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
needleDurations[1] = {"p_value": 1.0, "min": 1, "max": 6}
needleDurations[2] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[3] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[4] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
needleDurations[5] = {"min": 1, "max": 6}

PartnershipDurations = {"SEX": sexualDurations, "NEEDLE": needleDurations}

"""
Partnership acts distribution bins
"""
# REVIEW: Is this current logic?
sexualFrequency = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
sexualFrequency[1] = {"p_value": (0.323 + 0.262), "min": 1, "max": 6}
sexualFrequency[2] = {"p_value": (0.323 + 0.262 + 0.116), "min": 7, "max": 12}
sexualFrequency[3] = {"p_value": (0.323 + 0.262 + 0.116 + 0.121), "min": 13, "max": 24}
sexualFrequency[4] = {
    "p_value": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
sexualFrequency[5] = {
    "p_value": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
sexualFrequency[6] = {
    "p_value": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
sexualFrequency[7] = {
    "p_value": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
sexualFrequency[8] = {
    "p_value": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
sexualFrequency[9] = {"min": 37, "max": 48}

needleFrequency = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
needleFrequency[1] = {"p_value": 1.0, "min": 1, "max": 6}
needleFrequency[2] = {"p_value": (0.323 + 0.262), "min": 3, "max": 12}
needleFrequency[3] = {"p_value": (0.323 + 0.262), "min": 6, "max": 24}
needleFrequency[4] = {"p_value": (0.323 + 0.262), "min": 9, "max": 36}
needleFrequency[5] = {"min": 12, "max": 60}

PartnershipFrequency = {"SEX": sexualFrequency, "NEEDLE": needleFrequency}


"""
Sexual and injection transmission probabilities
"""
# REVIEW: is this the current logic?
SexTrans = {"MSM": {}, "HM": {}, "HF": {}}
SexTrans["MSM"] = {
    "0": 0.00745,
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

NeedleTrans = {
    "0": 0.007,
    "1": 0.007,
    "2": 0.0056,
    "3": 0.0028,
    "4": 0.0014,
    "5": 0.0002,
}
TransmissionProbabilities = {"SEX": SexTrans, "NEEDLE": NeedleTrans}


"""
Age bin distributions and HIV if utilized
"""
ageMatrix = {"WHITE": {}, "BLACK": {}}
ageMatrix["WHITE"] = {
    "Prop": {
        0: 0.0,
        1: 0.18,
        2: 0.18 + 0.16,
        3: 0.18 + 0.16 + 0.15,
        4: 0.18 + 0.16 + 0.15 + 0.20,
        5: 0.18 + 0.16 + 0.15 + 0.20 + 0.31,
    }
}
ageMatrix["BLACK"] = {
    "Prop": {
        0: 0.0,
        1: 0.28,
        2: 0.28 + 0.24,
        3: 0.28 + 0.24 + 0.19,
        4: 0.28 + 0.24 + 0.19 + 0.15,
        5: 0.28 + 0.24 + 0.19 + 0.15 + 0.14,
    }
}


"""
Clinic bins for targetting strategies
Bins represent partner numbers of the following category 0:0, 1:1, 2:2,  3:3-4, 4:5-9, 5:10+
"""
# REVIEW: clinic model might be completely out-of-date
clinicAgents = {"Low": {}, "Mid": {}, "High": {}}
clinicAgents["Mid"] = {
    0: {"Prob": 0.0, "min": 0, "max": 0},
    1: {"Prob": 0.054, "min": 0, "max": 1},
    2: {"Prob": 0.061, "min": 1, "max": 2},
    3: {"Prob": 0.168, "min": 3, "max": 4},
    4: {"Prob": 0.246, "min": 5, "max": 9},
    5: {"Prob": 0.471, "min": 10, "max": 120},
}
