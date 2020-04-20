from dotmap import DotMap

params = DotMap()
__author__ = "MaximilianKing"

from typing import Sequence, List, Dict, Optional, Any
from copy import deepcopy

"""
Main model parameters.
"""

####################
params.model.processes = 1  # number of processes in parallel (quadcore)
rSeed = (
    0  # seed for random number generator (0 for pure random, -1 for stepwise up to N_NC
)
params.model.seed.ppl = 0
params.model.seed.net = 0
params.model.seed.run = 0
N_MC = 1  # total number of iterations (Monte Carlo runs)
N_REPS = 1
params.model.num_pop = 5578  # population size
params.model.time_range = 24  # total time steps to iterate
params.model.burn_duration = 6
model = "Custom"  # Model Type for fast flag toggling
setting = "AtlantaMSM"
params.model.network.type = "max_k_comp_size"
####################

"""
Output flags and settings
"""
outputDir = ""

startAgentList = False
endingAgentList = False
intermAgentList = False
params.outputs.print_frequency = 12
MSMreport = True
HMreport = False
HFreport = False
params.outputs.draw_figures = False
drawEdgeList = True
params.outputs.calc_component_stats = True
params.features.agent_zero = False

params.outputs.reports = [
    "prepReport",
    "basicReport",
]

"""
Calibration scaling parameters for fitting to empirical data
"""

PARTNERTURNOVER = 0.9  # Partner acquisition parameters (higher number more partnering)
params.calibration.needle.sharing = 1.0  # PWID partner number scaling
params.calibration.needle.act = 1.0  # PWID act frequency scaling factor
params.calibration.sex.partner = 1.0  # Sexual partner number scaling factor
params.calibration.sex.act = 1.0  # Sexual acts  scaling factor
params.calibration.transmission = (
    1.0  # 0.92 # Global transmission probability scaling factor
)
params.hiv.acute.infectivity = (
    1.0  # Infectivity multiplier ratio for Acute status infections
)
params.hiv.dx.risk_reduction = (
    1.0  # Risk reduction in transmission probability for agents diagnosed
)
params.haart.transmission.prob = (
    1.0  # Scaling factor for effectiveness of ART therapy on xmission P
)
params.calibration.test_frequency = 1.0  # Scaling factor for testing frequency
params.calibration.mortality = 1.0  # Scaling factor for all cause mortality rates
params.hiv.aids.prob = 1.0  # Scaling factor for all progression to AIDS from HIV rates
params.calibration.haart_coverage = (
    1.0  # Scaling factor for enrollment on ART probability
)
params.calibration.incarceration = (
    1.0  # Scaling factor for probability of becoming incarcerated
)
params.demographics[
    "BLACK"
].transmission = (
    1.0  # Scaling factor for increased STI transmission P comparing race1/race2
)
params.calibration.partner_sample_depth = (
    100  # Sampling depth for partnering algorithm.
)

"""
Bond Params
"""
params.classes.bond_types_probs = {
    "social": 0.308,
    "multiplex": 0.105,
    "sexOnly": 0.587,
}
params.classes.bond_types_probs_PWID = {
    "social": 0.308,
    "multiplex": 0.105,
    "sexOnly": 0.587,
}
params.classes.bond_types = ["social"]
params.model.population.num_partners.type = "bins"

"""
Peer change params
"""
attitude = {0: 0.167, 1: 0.082, 2: 0.184, 3: 0.139, 4: 0.429}
PCA_PrEP = 0.30 * (
    1 - 0.575
)  # chance of attempting * chance of initiating oral or inj PrEP
opinion_threshold = (
    3.0  # opinion needed to initiate PrEP on a 0-4 scale (translated from 1-5 scale)
)
pcaChoice = "eigenvector"  # eigenvector or bridge, how the PCA is selected
awarenessProb = 0.055  # static probability of becoming spontaneously aware of PrEP
starting_awareness = 0.00  # awareness of PrEP at t0
knowledgeTransmission = (
    0.01  # per-act probability of knowledge change in unaware partner
)
opinionTransmission = (
    0.005  # per-act probability of opinion change in less-prominent partner
)
interactionProb: Dict[str, Any] = {
    "sexOnly": {1: {}},
    "multiplex": {1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}},
    "social": {1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}},
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


"""
High risk params
"""
params.high_risk.partner_scale = (
    300  # Linear increase to partner number during HR period
)
params.high_risk.proportion = (
    0.3  # Proportion of people who enter HR group when partner incarcerated
)
params.high_risk.sex_based.HM.duration = 6  # Duration of high risk for males
params.high_risk.sex_based.HF.duration = 6  # Duration of high risk for females
params.high_risk.condom_use_type = "Race"  # Race or Acts
HIV_MSMW = 0.0

# Misc. params
params.features.assort_mix = True
params.assort_mix.type = "Race"
flag_RaceAssortMix = True
params.assort_mix.coefficient = (
    0.75  # Proportion of race1 mixing with race2 when partnering.
)
params.needle_exchange.prevalence = 1.0  # Prevalence scalar on SSP
params.needle_exchange.init_at_pop = 0
params.needle_exchange.coverage = 0.0
params.model.network.component_size.max = 100

"""
Vaccine params
"""
params.vaccine.type = "RV144"
booster = False
params.vaccine.start = 1

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
inc_PtnrDissolution = 0.55
inc_treatment_startdate = 48  # Timestep where inc treatment can begin
inc_treatment_dur = (
    12  # Duration for which agents are forced on respective treatment post release
)
inc_treat_set = ["HM"]  # Set of agent classifiers effected by HR treatment
inc_treat_HRsex_beh = True  # Remove sexual higrisk behaviour during treatment duration
inc_treat_PWID_beh = True  # Remove PWID behav:iour during treatment duration
inc_treat_RIC = False  # Force retention in care of ART therapy

# PrEP params
init_with_vaccine = True
params.prep.type = ["Oral", "Inj"]  # Oral/Inj PrEP and/or vaccine
LAI_chance = 0.5
params.prep.target = (
    0.088  # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
)
params.prep.start = 0  # Start date for PrEP program (0 for start of model)
params.PWID.prep.adherence = 0.82  # Probability of being adherent
params.prep.efficacy.adherent = 0.96  # Efficacy of adherence PrEP
params.prep.efficacy.non_adherant = 0.76  # Efficacy of non-adherence PrEP
PrEP_falloutT = 0  # During PrEP remains effective post discontinuation
PrEP_resist = 0.01
params.prep.discontinue = 0.15
params.prep.target_model = "RandomTrial"  # Allcomers, Clinical
PrEP_clinic_cat = "Mid"  # If clinical target model, which category does it follow

if "Oral" in params.prep.type:
    params.PWID.prep.adherence = 1.0
    params.prep.efficacy.adherent = 0.96
    params.prep.efficacy.non_adherant = 0.76
    PrEP_falloutT = 1
    params.prep.discontinue = 0.00  # 0.15
    params.prep.peak_load = 1.0
    params.prep.half_life = 1.0
elif (
    "Inj" in params.prep.type
):  # TODO make both of these compatible; can use both at once??
    params.PWID.prep.adherence = 1.0
    params.prep.efficacy.adherent = 1.0
    params.prep.efficacy.non_adherant = 1.00
    PrEP_falloutT = 12
    params.prep.discontinue = 0.00  # 0.04
    params.prep.peak_load = 4.91
    params.prep.half_life = 40.0

"""
Model Type for fast flag toggling
    params.features.incar      Incarceration effects
    params.features.prep       PrEP enrollment
    flag_high_risk         High risk behavior for incar or genPop
    params.features.haart        ART therapy enrollment
    params.features.die_and_replace      Die and replace functionality

"""

####################

####################

if model == "PrEP":
    params.features.incar = False
    params.features.prep = True
    flag_high_risk = False
    params.features.haart = True
    params.features.die_and_replace = True
    params.features.static_network = False
    params.features.pca = True
elif model == "Incar":
    params.features.incar = True
    params.features.prep = False
    flag_high_risk = True
    params.features.haart = True
    params.features.die_and_replace = True
    params.features.static_network = False
    params.features.pca = True
elif model == "NoIncar":
    params.features.incar = False
    params.features.prep = False
    flag_high_risk = False
    params.features.haart = True
    params.features.die_and_replace = True
    params.features.static_network = False
    params.features.pca = True
elif model == "VaccinePrEP":
    params.features.incar = False
    params.features.prep = True
    flag_high_risk = False
    params.features.haart = True
    params.features.die_and_replace = True
    params.features.static_network = False
    params.vaccine.booster = True
    params.features.pca = True
elif model == "Custom":
    params.features.incar = False
    params.features.prep = True
    flag_high_risk = False
    params.features.haart = False
    params.features.die_and_replace = False
    params.features.static_network = True
    params.vaccine.booster = False
    params.features.pca = True

params.classes.sex_types = ["MSM"]
agentPopulations = deepcopy(params.classes.sex_types)
agentPopulations.append("PWID")

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
    "SAFESEX": 0.0,  # Probability of engaging in unsafe sex (per act)
    "NEEDLESH": 0.0,  # Probability of sharing syringes during join drug use (per act)
    "HIVTEST": 0.0,  # Probability of testing for HIV
    "INCAR": 0.0,  # Probability of becoming incarcerated (rate)
    "HAARTprev": 0.0,
    "HAARTadh": 0.0,  # Adherence to ART therapy
    "HAARTdisc": 0.0,  # Probability of discontinuing ART therapy
    "EligSE_PartnerType": [],  # List of agent SO types the agent cant partner with
    "PrEPdisc": 0.0,  # Probability of discontinuing PrEP treatment
    "HighRiskPrev": 0.0,
    "EligSE_PartnerType": [],
    "PrEPadh": 1.0,
    "PrEP_coverage": 0,
    "boosterInterval": 0,
    "boosterProb": 0,
    "vaccinePrev": 0,
    "vaccineInit": 0,
}

RaceClass1: Dict[str, Any] = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
RaceClass2: Dict[str, Any] = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
for a in ["MSM", "HM", "HF", "PWID"]:
    RaceClass1[a] = dict(RC_template)
    RaceClass2[a] = dict(RC_template)

RaceClass1["MSM"]["POP"] = 1.0
RaceClass1["MSM"]["HIV"] = 0.4
# StratW['MSM'] = {'POP':0.035, 'HIV':0.132, 'AIDS':0.048, 'HAARTprev':0.57, 'INCARprev':0.005, 'TestedPrev':0.84}

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
    {
        "Proportion": 0.00,
        "HAARTdisc": 0.018,
        "PrEPdisc": 0.0,
        "params.assort_mix.coefficient": 0.722,
    }
)

# RaceClass2 = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}, 'ALL':{}}
RaceClass2["MSM"].update(
    {
        "POP": 1.00,  # 0.028,
        "HIV": 0.367,
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
        "PrEPdisc": 0.00,
        "EligSE_PartnerType": "MSM",
        "PrEPadh": 0.568,
        "PrEP_coverage": 0.5,
        "vaccinePrev": 1.0,
        "boosterInterval": 3,
        "boosterProb": 1.0,
    }
)

RaceClass2["ALL"].update(
    {
        "Proportion": 1.00,
        "HAARTdisc": 0.018,
        "PrEPdisc": 0.0,
        "params.assort_mix.coefficient": 0.765,
    }
)

params.demographics = {"WHITE": RaceClass1, "BLACK": RaceClass2}

"""
Partnership durations and
"""
params.partnership.sex.duration: Dict[int, Any] = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
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

params.partnership.needle.duration: Dict[int, Any] = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
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
Partnership acts and
"""
params.partnership.sex.frequency: Dict[int, Any] = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
params.partnership.sex.frequency[1] = {"prob": (0.323 + 0.262), "min": 1, "max": 6}
params.partnership.sex.frequency[2] = {
    "prob": (0.323 + 0.262 + 0.116),
    "min": 7,
    "max": 12,
}
params.partnership.sex.frequency[3] = {
    "prob": (0.323 + 0.262 + 0.116 + 0.121),
    "min": 13,
    "max": 24,
}
params.partnership.sex.frequency[4] = {
    "prob": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
params.partnership.sex.frequency[5] = {
    "prob": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
params.partnership.sex.frequency[6] = {
    "prob": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
params.partnership.sex.frequency[7] = {
    "prob": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
params.partnership.sex.frequency[8] = {
    "prob": (0.323 + 0.262 + 0.116 + 0.121 + 0.06),
    "min": 25,
    "max": 36,
}
params.partnership.sex.frequency[9] = {"min": 37, "max": 48}

params.partnership.needle.frequency: Dict[int, Any] = {
    1: {},
    2: {},
    3: {},
    4: {},
    5: {},
}
params.partnership.needle.frequency[1] = {"prob": 1.0, "min": 1, "max": 6}
params.partnership.needle.frequency[2] = {"prob": (0.323 + 0.262), "min": 1, "max": 6}
params.partnership.needle.frequency[3] = {"prob": (0.323 + 0.262), "min": 1, "max": 6}
params.partnership.needle.frequency[4] = {"prob": (0.323 + 0.262), "min": 1, "max": 6}
params.partnership.needle.frequency[5] = {"min": 1, "max": 6}

PartnershipDurations = {
    "SEX": params.partnership.sex.duration,
    "NEEDLE": params.partnership.needle.duration,
}

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

params.partnership.needle.transmission = {
    "0": 0.007,
    "1": 0.007,
    "2": 0.0056,
    "3": 0.0028,
    "4": 0.0014,
    "5": 0.0002,
}

TransmissionProbabilities: Dict[str, Any] = {
    "SEX": SexTrans,
    "NEEDLE": params.partnership.needle.transmission,
}

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
