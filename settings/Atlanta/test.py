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
N_POP = 17440  # population size
TIME_RANGE = 120  # total time steps to iterate
burnDuration = 36
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
drawEdgeList = False

reports = [
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
cal_ProgAIDS = 0.05  # Scaling factor for all progression to AIDS from HIV rates
cal_ART_cov = 0.2  # Scaling factor for enrollment on ART probability
cal_IncarP = 1.0
cal_raceXmission = 1.95
cal_ptnrSampleDepth = 100
cal_Vaccine = 0  # determines vaccine initiation during run

"""
Network Params
"""
bond_type = []
mean_partner_type = "mean"

"""
Peer Change Params
"""
flag_PCA = False

"""
High risk params
"""
HR_partnerScale = 300  # Linear increase to partner number during HR period
HR_proportion = 0.3  # Proportion of people who enter HR group when partner incarcerated
HR_M_dur = 6  # Duration of high risk for males
HR_F_dur = 6  # Duration of high risk for females
condomUseType = "Race"  # Race or Acts

"""
Misc. params
"""
flag_AssortativeMix = True
AssortMixType = "Race"
flag_RaceAssortMix = True
AssortMixCoeff = 0.75  # Proportion of race1 mixing with race2 when partnering.
safeNeedleExchangePrev = 1.0  # Prevalence scalar on SNE
initTreatment = 0
treatmentCov = 0.0

"""
Vaccine params
"""
vaccine_type = "HVTN702"
booster = False
vaccine_start = 1
init_with_vaccine = False

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
PrEP_type = ["Oral", "Vaccine"]  # Oral/Inj PrEP and/or vaccine
PrEP_Target = (
    1.0  # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
)
PrEP_startT = 0  # Start date for PrEP program (0 for start of model)
PrEP_Adherence = 0.82  # Probability of being adherent
PrEP_AdhEffic = 0.96  # Efficacy of adherence PrEP
PrEP_NonAdhEffic = 0.76  # Efficacy of non-adherence PrEP
PrEP_falloutT = 0  # During PrEP remains effective post discontinuation
PrEP_resist = 0.01
PrEP_disc = 0.15
PrEP_target_model = {"Racial"}
PrEP_init_var1 = 0.5
PrEP_init_var2 = 0.05
PrEP_clinic_cat = ""

if "Oral" in PrEP_type:
    PrEP_Adherence = "AtlantaMSM"
    PrEP_AdhEffic = 0.96
    PrEP_NonAdhEffic = 0.76
    PrEP_falloutT = 1
    PrEP_disc = 0.15
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
    flag_high_risk         High risk behavior for incar or genPop
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
elif model == "Incar":
    flag_incar = True
    flag_PrEP = False
    flag_high_risk = True
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
elif model == "NoIncar":
    flag_incar = False
    flag_PrEP = False
    flag_high_risk = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
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

agentPopulations = ["MSM", "HF", "HM", "IDU"]
agentSexTypes = ["MSM", "HF", "HM", "IDU"]

"""
RaceClass is a distinct racial/ethnic/social classification for demographics of the population.
"""