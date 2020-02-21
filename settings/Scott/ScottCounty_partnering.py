from dotmap import DotMap

params = DotMap()
__author__ = "MaximilianKing"


"""
Main model parameters.
"""

####################
params.model.processes = 1  # number of processes in parallel (quadcore)
params.model.seed.ppl = (
    0  # seed for random number generator (0 for pure random, -1 for stepwise up to N_NC
)
params.model.seed.net = 0
params.model.seed.run = 0
N_MC = 1  # total number of iterations (Monte Carlo runs)
N_REPS = 1
params.model.num_pop = 1000  # population size
params.model.time_range = 120  # total time steps to iterate
params.model.burn_duration = 36
model = "StaticZero"  # Model Type for fast flag toggling
params.model.network.type = "scale_free"  # scale_free or max_k_comp_size
setting = "Scott"
####################


"""
Output flags and settings
"""
outputDir = ""

printIncidenceEvents = False
printStartAgentList = False
printEndingAgentList = False
printIntermAgentList = False
params.outputs.print_frequency = 60
params.outputs.calc_network_stats = False
params.outputs.calc_component_stats = False
params.outputs.draw_figures = False
drawEdgeList = True
drawFigureColor = "MSW"

params.outputs.reports = [
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
params.calibration.needle.sharing = 1.0  # PWID partner number scaling
params.calibration.needle.act = 2.0  # PWID act frequency scaling factor
params.calibration.sex.partner = 1.0  # Sexual partner number scaling factor
params.calibration.sex.act = 1.0  # Sexual acts  scaling factor
params.calibration.transmission = 1.0  # Global transmission probability scaling factor
params.calibration.acute = (
    4.3  # Infectivity multiplier ratio for Acute status infections
)
params.calibration.risk_reduction.transmission = (
    0.53  # Risk reduction in transmission probability for agents diagnosed
)
params.calibration.risk_reduction.haart = (
    1.0  # Scaling factor for effectiveness of ART therapy on xmission P
)
params.calibration.test_frequency = 1.0  # Scaling factor for testing frequency
params.calibration.mortality = 0.5  # Scaling factor for all cause mortality rates
params.calibration.aids_progression = (
    1.0  # Scaling factor for all progression to AIDS from HIV rates
)
params.calibration.haart_coverage = (
    1.0  # Scaling factor for enrollment on ART probability
)
params.calibration.incarceration = 1.0
params.calibration.race_transmission = 1.0
params.calibration.partner_sample_depth = 100

"""
Bond Params
"""
params.classes.bond_types = []
params.model.population.num_partners.type = "mean"

"""
Peer Change Params
"""
params.features.pca = False


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
params.high_risk.condom_use_type = "Acts"  # race or number of acts

"""
Misc. params
"""

params.features.assort_mix = False  # Boolean for if assortative mixing occurs at all
params.assort_mix.type = None  # Other assortative mixing types
flag_RaceAssortMix = False  # Assortative mix by race
params.assort_mix.coefficient = 0.8  # Proportion of following given assort mix rules
params.needle_exchange.prevalence = 1.0  # Prevalence scalar on SNE
params.needle_exchange.init_at_pop = 10  # Requirement to start treatment
params.needle_exchange.coverage = 0.60  # Prop that receive treatment
limitComponentSize = False
params.model.network.component_size.max = 100
params.model.network.component_size.min = 2

"""
Incarceration params
"""
inc_JailMax = 22
inc_JailMin = 8
inc_JailTestProb = 0.69
inc_PrisMax = 96
inc_PrisMin = 45
inc_PrisTestProb = 0.69
inc_PropPrison = 0.5
params.incar.haart.prob = 0.51
params.incar.haart.adherence = 0.21
params.incar.haart.discontinue = 0.12
inc_Recidivism = 0.267
inc_PtnrDissolution = 0.55
inc_treatment_dur = (
    6  # Duration for which agents are forced on respective treatment post release
)
inc_treat_set = ["HM"]  # Set of agent classifiers effected by HR treatment
inc_treat_behavior = True  # Remove PWID behaviour during treatment duration
inc_treat_RIC = False  # Force retention in care of ART therapy

"""
PrEP params
"""
params.prep.type = "Oral"  # Oral/Inj PrEP modes
params.prep.target = (
    0.000  # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
)
params.prep.start = 0  # Start date for PrEP program (0 for start of model)
params.PWID.prep.adherence = 0.82  # Probability of being adherent
params.prep.efficacy.adherent = 0.96  # Efficacy of adherence PrEP
params.prep.efficacy.non_adherant = 0.76  # Efficacy of non-adherence PrEP
PrEP_falloutT = 0  # During PrEP remains effective post discontinuation
PrEP_resist = 0.01
params.prep.discontinue = 0.15
params.prep.target_model = "Allcomers"  # Clinical, Allcomers
PrEP_clinic_cat = "Mid"

if params.prep.type == "Oral":
    params.PWID.prep.adherence = 0.923
    params.prep.efficacy.adherent = 0.96
    params.prep.efficacy.non_adherant = 0.76
    PrEP_falloutT = 1
    params.prep.discontinue = 0.15
elif params.prep.type == "Inj":
    params.PWID.prep.adherence = 1.0
    params.prep.efficacy.adherent = 1.0
    params.prep.efficacy.non_adherant = 1.00
    PrEP_falloutT = 12
    params.prep.discontinue = 0.04
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
if model == "PrEP":
    params.features.incar = False
    params.features.prep = True
    flag_high_risk = False
    params.features.haart = True
    params.features.die_and_replace = True
    params.features.static_n = False
    params.features.agent_zero = False

elif model == "Incar":
    params.features.incar = True
    params.features.prep = False
    flag_high_risk = True
    params.features.haart = True
    params.features.die_and_replace = True
    params.features.static_n = False
    params.features.agent_zero = False

elif model == "NoIncar":
    params.features.incar = False
    params.features.prep = False
    flag_high_risk = True
    params.features.haart = True
    params.features.die_and_replace = True
    params.features.static_n = False
    params.features.agent_zero = False

elif model == "StaticZero":
    params.features.incar = False
    params.features.prep = False
    flag_high_risk = False
    params.features.haart = False
    params.features.die_and_replace = False
    params.features.static_n = True
    params.features.agent_zero = True

elif model == "Custom":
    params.features.incar = False
    params.features.prep = False
    flag_high_risk = False
    params.features.haart = False
    params.features.die_and_replace = False
    params.features.static_n = True
    params.features.agent_zero = False

agentPopulations = ["HM", "HF"]
params.classes.sex_types = ["HM", "HF", "MSM"]
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
    "SAFESEX": 0.0,  # Probability of engaging in safe sex (per act)
    "NEEDLESH": 0.0,  # Probability of sharing syringes during join drug use (per act)
    "HIVTEST": 0.0,  # Probability of testing for HIV
    "INCAR": 0.0,  # Probability of becoming incarcerated (rate)
    "HAARTadh": 0.0,  # Adherence to ART therapy
    "HAARTdisc": 0.0,  # Probability of discontinuing ART therapy
    "PrEPprev": 0.0,  # Proportion of HIV- that are enrolled on PrEP
    "PrEPdisc": 0.0,  # Probability of discontinuing PrEP treatment
    "EligSE_PartnerType": None,  # List of agent SO types the agent cant partner with
    "AssortMixMatrix": [],  # List of assortMix Matrix to be zipped with EligPart
}

RC_allTemplate = {
    "Proportion": 1.00,  # Proportion of total population that is raceclass
    "HAARTdisc": 0.018,  # Overall HAART discontinuation probability
    "PrEPdisc": 0.0,  # Overall PrEP discontinuation probability
    "params.assort_mix.coefficient": 1.0,  # Proportion RC mixes with other raceclass
}

RaceClass1 = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
RaceClass2 = {"MSM": {}, "HM": {}, "HF": {}, "PWID": {}, "ALL": {}}
for a in ["MSM", "HM", "HF", "PWID"]:
    RaceClass1[a] = dict(RC_template)
    RaceClass2[a] = dict(RC_template)

RaceClass1["HM"].update(
    {
        "POP": 0.577,
        "HIV": 0.0014,
        "AIDS": 0.6780,
        "HAARTprev": 0.41,
        "INCARprev": 0.0274,
        "TestedPrev": 0.90,
        "HighRiskPrev": 0.0,
        "NUMPartn": 1.2,
        "NUMSexActs": 13.4,
        "SAFESEX": 0.89,
        "NEEDLESH": 0.43,
        "HIVTEST": 0.034,
        "INCAR": 0.00,
        "HAARTadh": 0.405,
        "HAARTdisc": 0.000,
        "PrEPdisc": 0.0000,
        "EligSE_PartnerType": "HF",
    }
)

RaceClass1["HF"].update(
    {
        "POP": 0.423,
        "HIV": 0.0004,
        "AIDS": 0.573,
        "HAARTprev": 0.47,
        "INCARprev": 0.000,
        "TestedPrev": 0.90,
        "HighRiskPrev": 0.0,
        "NUMPartn": 1.2,
        "NUMSexActs": 12.74,
        "SAFESEX": 0.43,
        "NEEDLESH": 0.43,
        "HIVTEST": 0.034,
        "INCAR": 0.00,
        "HAARTadh": 0.405,
        "HAARTdisc": 0.000,
        "PrEPdisc": params.prep.discontinue,
        "EligSE_PartnerType": "HM",
    }
)

RaceClass1["PWID"].update(
    {
        "POP": 0.847,
        "HIV": 0.000,
        "AIDS": 0.6780,
        "HAARTprev": 0.41,
        "INCARprev": 0.0274,
        "TestedPrev": 0.90,
        "HighRiskPrev": 0.0,
        "NUMPartn": 400.606,
        "NUMSexActs": 5.0,
        "SAFESEX": 0.89,
        "NEEDLESH": 0.63,
        "HIVTEST": 0.055,
        "INCAR": 0.00,
        "HAARTadh": 0.405,
        "HAARTdisc": 0.000,
        "PrEPadh": 0.55,
        "PrEPdisc": 0.0000,
        "EligSE_PartnerType": "PWID",
    }
)

RaceClass1["ALL"].update(
    {
        "Proportion": 1.00,
        "HAARTdisc": 0.018,
        "PrEPdisc": 0.0,
        "params.assort_mix.coefficient": 1.0,
    }
)

RaceClass2["ALL"].update(
    {
        "Proportion": 0.00,
        "HAARTdisc": 0.018,
        "PrEPdisc": 0.0,
        "params.assort_mix.coefficient": 1.0,
    }
)

params.demographics = {"WHITE": RaceClass1, "BLACK": RaceClass2}


"""
Partnership duration distribution bins
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
Partnership acts distribution bins
"""
params.partnership.sex.frequency = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
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

params.partnership.needle.frequency = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
params.partnership.needle.frequency[1] = {"prob": 1.0, "min": 1, "max": 6}
params.partnership.needle.frequency[2] = {"prob": (0.323 + 0.262), "min": 3, "max": 12}
params.partnership.needle.frequency[3] = {"prob": (0.323 + 0.262), "min": 6, "max": 24}
params.partnership.needle.frequency[4] = {"prob": (0.323 + 0.262), "min": 9, "max": 36}
params.partnership.needle.frequency[5] = {"min": 12, "max": 60}

PartnershipFrequency = {
    "SEX": params.partnership.sex.frequency,
    "NEEDLE": params.partnership.needle.frequency,
}


"""
Sexual and injection transmission probabilities
"""
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
    },
    "HIV": {0: 0.0, 1: 0.006, 2: 0.029, 3: 0.055, 4: 0.069, 5: 0.025},
}
ageMatrix["BLACK"] = {
    "Prop": {
        0: 0.0,
        1: 0.28,
        2: 0.28 + 0.24,
        3: 0.28 + 0.24 + 0.19,
        4: 0.28 + 0.24 + 0.19 + 0.15,
        5: 0.28 + 0.24 + 0.19 + 0.15 + 0.14,
    },
    "HIV": {0: 0.0, 1: 0.144, 2: 0.144, 3: 0.250, 4: 0.377, 5: 0.194},
}


"""
Age mixing matrix for assortative mixing by age
"""
mixingMatrix = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}}
mixingMatrix[1] = {1: 0.500, 2: 0.226, 3: 0.123, 4: 0.088, 5: 0.064}
mixingMatrix[2] = {1: 0.156, 2: 0.500, 3: 0.185, 4: 0.099, 5: 0.061}
mixingMatrix[3] = {1: 0.074, 2: 0.162, 3: 0.500, 4: 0.184, 5: 0.079}
mixingMatrix[4] = {1: 0.057, 2: 0.093, 3: 0.199, 4: 0.500, 5: 0.150}
mixingMatrix[5] = {1: 0.062, 2: 0.086, 3: 0.128, 4: 0.224, 5: 0.500}

"""
Clinic bins for targetting strategies
Bins represent partner numbers of the following category 0:0, 1:1, 2:2,  3:3-4, 4:5-9, 5:10+
"""
clinicAgents = {"Low": {}, "Mid": {}, "High": {}}
clinicAgents["Mid"] = {
    0: {"Prob": 0.0, "min": 0, "max": 0},
    1: {"Prob": 0.054, "min": 0, "max": 1},
    2: {"Prob": 0.061, "min": 1, "max": 2},
    3: {"Prob": 0.168, "min": 3, "max": 4},
    4: {"Prob": 0.246, "min": 5, "max": 9},
    5: {"Prob": 0.471, "min": 10, "max": 120},
}
