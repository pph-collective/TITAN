__author__ = 'MaximilianKing'


"""
Main model parameters.
"""

####################
PROCESSES = 1           # number of processes in parallel (quadcore)
rSeed = 0               # seed for random number generator (0 for pure random, -1 for stepwise up to N_NC
N_MC = 1              # total number of iterations (Monte Carlo runs)
N_POP = 1000          # population size
TIME_RANGE = 12        # total time steps to iterate
burnDuration = 3	# total time for burning in period (equillibration)
model = 'Incar'         # Model Type for fast flag toggling
setting = 'Cali'
####################

"""
Output flags and settings
"""
outputDir = ''

startAgentList = False
endingAgentList = False
intermAgentList = False
intermPrintFreq = 10
MSMreport = True
HMreport = False
HFreport = False
drawFigures = False


"""
Calibration scaling parameters for fitting to empirical data
"""

PARTNERTURNOVER = 0.2           # Partner acquisition parameters (higher number more partnering)
cal_NeedlePartScaling = 1.0     # IDU partner number scaling
cal_NeedleActScaling = 1.0      # IDU act frequency scaling factor
cal_SexualPartScaling = 1.0     # Sexual partner number scaling factor
cal_SexualActScaling = 1.0      # Sexual acts  scaling factor
cal_pXmissionScaling = 1.0      # Global transmission probability scaling factor
cal_AcuteScaling = 10.0         # Infectivity multiplier ratio for Acute status infections
cal_RR_Dx = 0.50                # Risk reduction in transmission probability for agents diagnosed
cal_RR_HAART = 1.0              # Scaling factor for effectiveness of ART therapy on xmission P
cal_TestFreq = 1.0              # Scaling factor for testing frequency
cal_Mortality = 0.5             # Scaling factor for all cause mortality rates
cal_ProgAIDS = 1.0              # Scaling factor for all progression to AIDS from HIV rates
cal_ART_cov = 1.0               # Scaling factor for enrollment on ART probability
cal_IncarP = 1.0                # Scaling factor for probability of becoming incarcerated
cal_raceXmission = 1.0          # Scaling factor for increased STI transmission P comparing race1/race2
cal_ptnrSampleDepth = 100       # Sampling depth for partnering algorithm.

"""
High risk params
"""
HR_partnerScale = 300           # Linear increase to partner number during HR period
HR_proportion = 0.3             # Proportion of people who enter HR group when partner incarcerated
HR_M_dur = 6                    # Duration of high risk for males
HR_F_dur = 6                    # Duration of high risk for females


"""
Misc. params
"""
flag_AssortativeMix = True     # Boolean for if assortative mixing occurs at all
AssortMixType = "HR"            # Other assortative mixing types
flag_AgeAssortMix = False       # Assortative mix by age
flag_RaceAssortMix = False      # Assortative mix by race
AssortMixCoeff = 0.3           # Proportion of following given assort mix rules
safeNeedleExchangePrev = 1.0    # Prevalence scalar on SNE

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
inc_ARTenroll = 0.51
inc_ARTadh = 0.21
inc_ARTdisc = 0.12
inc_Recidivism = 0.267
inc_PtnrDissolution = 0.55
inc_treatment_dur = 6           # Duration for which agents are forced onto treatment post release
inc_treat_set = ['HM','IDU']         # Set of agent classifiers effected by HR treatment
inc_treat_IDU = True

"""
PrEP params
"""
PrEP_type = "Oral"              #Oral/Inj PrEP modes
PrEP_Target = 10.4              # Target coverage for PrEP therapy at 10 years (unused in non-PrEP models)
PrEP_startT = 0                 # Start date for PrEP program (0 for start of model)
PrEP_Adherence = 0.82           # Probability of being adherent
PrEP_AdhEffic = 0.96            # Efficacy of adherence PrEP
PrEP_NonAdhEffic = 0.76         # Efficacy of non-adherence PrEP
PrEP_falloutT = 0               # During PrEP remains effective post discontinuation
PrEP_resist = 0.01              # Probability of PrEP resistance developing
PrEP_disc = 0.15                # Per month probability of PrEP discontinuation
PrEP_target_model = 'MSM'       # Allcomers, Clinical, Allcomers, HighPN5, HighPN10, SRIns, SR,Rec, MSM
PrEP_clinic_cat = 'Mid'         # If clinical target model, which category does it follow

if PrEP_type == 'Oral':
    PrEP_Adherence = 0.923
    PrEP_AdhEffic = 0.96
    PrEP_NonAdhEffic = 0.76
    PrEP_falloutT = 1
    PrEP_disc = 0.15
elif PrEP_type == 'Inj':
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
if model == 'PrEP':
    flag_incar = False
    flag_PrEP = True
    flag_HR = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_agentZero = False

elif model == 'Incar':
    flag_incar = True
    flag_PrEP = False
    flag_HR = True
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_agentZero = False

elif model == 'NoIncar':
    flag_incar = False
    flag_PrEP = False
    flag_HR = True
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_agentZero = False

elif model == 'StaticZero':
    flag_incar = False
    flag_PrEP = False
    flag_HR = False
    flag_ART = False
    flag_DandR = False
    flag_staticN = True
    flag_agentZero = False

elif model == 'Custom':
    flag_incar = False
    flag_PrEP = True
    flag_HR = False
    flag_ART = True
    flag_DandR = True
    flag_staticN = False
    flag_agentZero = False

agentSexTypes = ['HM', 'HF', 'MSM', 'MTF']
"""
RaceClass is a distinct racial/ethnic/social classification for demographics of the population.
ID is the specific mode of partnership the agent engages in (ie MSM, HM, HF, PWID)
RaceClass agent classifier template
"""
RC_template = {     'Race':None,            #Race of demographic
                    'Class':None,           #Classification of networking
                    'POP':0.0,              #Percentage of total raceclass population that are ID
                    'HIV':0.0,              #Proportion of total ID population that are HIV
                    'AIDS':0.0,             #Proportion of total HIV_ID that are AIDS
                    'HAARTprev':0.0,        #Proportion of HIV_TESTED_ID that are enrolled on ART
                    'INCARprev':0.0,        #Proportion of ID that are incarcerated
                    'TestedPrev':0.0,       #Proportion of HIV_ID that are tested positively
                    'mNPart':0.0,           #Mean number of sex partners
                    'NUMPartn':0.0,         #Number of partners (redundant)
                    'NUMSexActs':0.0,       #Mean number of sex acts with each partner
                    'UNSAFESEX':0.0,        #Probability of engaging in unsafe sex (per act)
                    'NEEDLESH':0.0,         #Probability of sharing syringes during join drug use (per act)
                    'HIVTEST':0.0,          #Probability of testing for HIV
                    'INCAR':0.0,            #Probability of becoming incarcerated (rate)
                    'HAARTadh':0.0,         #Adherence to ART therapy
                    'HAARTdisc':0.0,        #Probability of discontinuing ART therapy
                    'PrEPdisc':0.0,         #Probability of discontinuing PrEP treatment
                    'EligPartnerType':[],   #List of agent SO types the agent cant partner with
                    'AssortMixMatrix':[]    #List of assortMix Matrix to be zipped with EligPart
                }

RC_allTemplate = {  'Proportion':1.00,      #Proportion of total population that is raceclass
                    'HAARTdisc':0.018,      #Overall HAART discontinuation probability
                    'PrEPdisc':0.0,         #Overall PrEP discontinuation probability
                    'AssortMixCoeff':1.0,   #Proportion RC mixes with other raceclass
                }

RaceClass1 = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}, 'ALL':{}}
RaceClass2 = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}, 'ALL':{}}
for a in ['MSM','HM','HF','PWID']:
    RaceClass1[a] = dict(RC_template)
    RaceClass2[a] = dict(RC_template)

RaceClass1['HM'] = {'POP':0.4224565212,
                     'HIV':0.154,
                     'AIDS':0.189,
                     'HAARTprev':0.919,
                     'INCARprev':0.100,
                     'TestedPrev':0.347,
                     'NUMPartn':3.0,
                     'NUMSexActs':3.4,
                     'UNSAFESEX':0.77,
                     'NEEDLESH':0.00,
                     'HIVTEST':0.035,
                     'INCAR':0.100,
                     'HAARTadh':0.67,
                     'HAARTdisc':0.000,
                     'PrEPadh':0.55,
                     'PrEPdisc':0.0000,
                     'EligPartnerType':['HF']
                     }

RaceClass1['HF'] = {'POP':0.4993509943,
                     'HIV':0.21,
                     'AIDS':0.205,
                     'HAARTprev':0.859,
                     'INCARprev':0.000,
                     'TestedPrev':0.653,
                     'NUMPartn':2.0,
                     'NUMSexActs':4.6,
                     'UNSAFESEX':0.67,
                     'NEEDLESH':0.00,
                     'HIVTEST':0.031,
                     'INCAR':0.00,
                     'HAARTadh':0.62,
                     'HAARTdisc':0.000,
                     'PrEPadh':0.55,
                     'PrEPdisc':PrEP_disc,
                     'EligPartnerType':['HM']
                     }

RaceClass1['MSM'] = {'POP':0.07689447306,
                     'HIV':0.2093,
                     'AIDS':0.079,
                     'HAARTprev':0.926,
                     'INCARprev':0.000,
                     'TestedPrev':0.956,
                     'NUMPartn':4.0,
                     'NUMSexActs':2.8,
                     'UNSAFESEX':0.49,
                     'NEEDLESH':0.00,
                     'HIVTEST':0.13,
                     'INCAR':0.00,
                     'HAARTadh':0.66,
                     'HAARTdisc':0.000,
                     'PrEPadh':0.55,
                     'PrEPdisc':PrEP_disc,
                     'EligPartnerType':['MSM']
                     }

RaceClass1['MTF'] = {'POP':0.00129801141,
                     'HIV':0.33986,
                     'AIDS':0.636,
                     'HAARTprev':1.00,
                     'INCARprev':0.000,
                     'TestedPrev':0.95,
                     'NUMPartn':4.7,
                     'NUMSexActs':4.6,
                     'UNSAFESEX':0.644,
                     'NEEDLESH':0.00,
                     'HIVTEST':0.155,
                     'INCAR':0.00,
                     'HAARTadh':0.67,
                     'HAARTdisc':0.000,
                     'PrEPadh':0.55,
                     'PrEPdisc':PrEP_disc,
                     'EligPartnerType':['MSM']
                     }

RaceClass1['PWID'] = {'POP':0.026,
                     'HIV':0.146,
                     'AIDS':0.1765,
                     'HAARTprev':0.609,
                     'INCARprev':0.00,
                     'TestedPrev':0.61,
                     'NUMPartn':2.0,
                     'NUMSexActs':3.5,
                     'UNSAFESEX':0.00,
                     'NEEDLESH':0.092,
                     'HIVTEST':0.077,
                     'INCAR':0.001,
                     'HAARTadh':0.605,
                     'HAARTdisc':0.000,
                     'PrEPadh':0.25,
                     'PrEPdisc':0.0000,
                     'EligPartnerType':['IDU']
                     }


RaceClass1['ALL'] = {'Proportion':1.00,
                      'HAARTdisc':0.018,
                     'PrEPdisc':0.0,
                     'AssortMixCoeff':1.0,
                      }

RaceClass2['ALL'] = {'Proportion':0.00,
                      'HAARTdisc':0.018,
                     'PrEPdisc':0.0,
                     'AssortMixCoeff':1.0,
                      }

DemographicParams = {'WHITE':RaceClass1, 'BLACK':RaceClass2}


"""
Partnership duration distribution bins
"""
sexualDurations = {1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
sexualDurations[1] = {'p_value':(0.323 + 0.262), 'min':1, 'max':6}
sexualDurations[2] = {'p_value':(0.323 + 0.262 + 0.116), 'min':7, 'max':12}
sexualDurations[3] = {'p_value':(0.323 + 0.262 + 0.116 + 0.121), 'min':13, 'max':24}
sexualDurations[4] = {'p_value':(0.323 + 0.262 + 0.116 + 0.121 + 0.06), 'min':25, 'max':36}
sexualDurations[5] = {'min':37, 'max':48}

needleDurations = {1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
needleDurations[1] = {'p_value':1.0, 'min':1, 'max':6}
needleDurations[2] = {'p_value':(0.323 + 0.262), 'min':1, 'max':6}
needleDurations[3] = {'p_value':(0.323 + 0.262), 'min':1, 'max':6}
needleDurations[4] = {'p_value':(0.323 + 0.262), 'min':1, 'max':6}
needleDurations[5] = {'min':1, 'max':6}

PartnershipDurations = {'SEX':sexualDurations, 'NEEDLE':needleDurations}

"""
Partnership acts distribution bins
"""
sexualFrequency = {1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
sexualFrequency[1] = {'p_value':(0.323 + 0.262), 'min':1, 'max':6}
sexualFrequency[2] = {'p_value':(0.323 + 0.262 + 0.116), 'min':7, 'max':12}
sexualFrequency[3] = {'p_value':(0.323 + 0.262 + 0.116 + 0.121), 'min':13, 'max':24}
sexualFrequency[4] = {'p_value':(0.323 + 0.262 + 0.116 + 0.121 + 0.06), 'min':25, 'max':36}
sexualFrequency[5] = {'p_value':(0.323 + 0.262 + 0.116 + 0.121 + 0.06), 'min':25, 'max':36}
sexualFrequency[6] = {'p_value':(0.323 + 0.262 + 0.116 + 0.121 + 0.06), 'min':25, 'max':36}
sexualFrequency[7] = {'p_value':(0.323 + 0.262 + 0.116 + 0.121 + 0.06), 'min':25, 'max':36}
sexualFrequency[8] = {'p_value':(0.323 + 0.262 + 0.116 + 0.121 + 0.06), 'min':25, 'max':36}
sexualFrequency[9] = {'min':37, 'max':48}

needleFrequency = {1:{}, 2:{}, 3:{}, 4:{}, 5:{}}
needleFrequency[1] = {'p_value':1.0, 'min':1, 'max':6}
needleFrequency[2] = {'p_value':(0.323 + 0.262), 'min':3, 'max':12}
needleFrequency[3] = {'p_value':(0.323 + 0.262), 'min':6, 'max':24}
needleFrequency[4] = {'p_value':(0.323 + 0.262), 'min':9, 'max':36}
needleFrequency[5] = {'min':12, 'max':60}

PartnershipFrequency = {'SEX':sexualFrequency, 'NEEDLE':needleFrequency}



"""
Sexual and injection transmission probabilities
"""
SexTrans = {'MSM':{}, 'HM':{}, 'HF':{}}
SexTrans['MSM'] = {'0':0.00745, '1':0.005, '2':0.004, '3':0.002, '4':0.001, '5':0.0001}
SexTrans['MTF'] = {'0':0.00745, '1':0.005, '2':0.004, '3':0.002, '4':0.001, '5':0.0001}
SexTrans['HM'] = {'0':0.001, '1':0.001, '2':0.0008, '3':0.0004, '4':0.0002, '5':0.0001}
SexTrans['HF'] = {'0':0.001, '1':0.001, '2':0.0008, '3':0.0004, '4':0.0002, '5':0.0001}

NeedleTrans = {'0':0.007, '1':0.007, '2':0.0056, '3':0.0028, '4':0.0014, '5':0.0002}
TransmissionProbabilities = {'SEX':SexTrans, 'NEEDLE':NeedleTrans}



"""
Age bin distributions and HIV if utilized
"""
ageMatrix = {'WHITE':{}, 'BLACK':{}}
ageMatrix['WHITE'] = {'Prop':{0:0.0, 1:0.18, 2:0.18+0.16, 3:0.18+0.16+0.15, 4:0.18+0.16+0.15+0.20, 5:0.18+0.16+0.15+0.20+0.31},
                       'HIV':{0:0.0, 1:0.006, 2:0.029, 3:0.055, 4:0.069, 5:0.025}}
ageMatrix['BLACK'] = {'Prop':{0:0.0, 1:0.28, 2:0.28+0.24, 3:0.28+0.24+0.19, 4:0.28+0.24+0.19+0.15, 5:0.28+0.24+0.19+0.15+0.14},
                       'HIV':{0:0.0, 1:0.144, 2:0.144, 3:0.250, 4:0.377, 5:0.194}}


"""
Age mixing matrix for assortative mixing by age
"""
mixingMatrix = {1:{},2:{},3:{},4:{},5:{}}
mixingMatrix[1] = {1:0.500,2:0.226,3:0.123,4:0.088,5:0.064}
mixingMatrix[2] = {1:0.156,2:0.500,3:0.185,4:0.099,5:0.061}
mixingMatrix[3] = {1:0.074,2:0.162,3:0.500,4:0.184,5:0.079}
mixingMatrix[4] = {1:0.057,2:0.093,3:0.199,4:0.500,5:0.150}
mixingMatrix[5] = {1:0.062,2:0.086,3:0.128,4:0.224,5:0.500}

"""
Clinic bins for targetting strategies
Bins represent partner numbers of the following category 0:0, 1:1, 2:2,  3:3-4, 4:5-9, 5:10+
"""
clinicAgents = {'Low':{},'Mid':{},'High':{}}
clinicAgents['Mid'] = {0:{'Prob':0.0,'min':0,'max':0},
                       1:{'Prob':0.054,'min':0,'max':1},
                       2:{'Prob':0.061,'min':1,'max':2},
                       3:{'Prob':0.168,'min':3,'max':4},
                       4:{'Prob':0.246,'min':5,'max':9},
                       5:{'Prob':0.471,'min':10,'max':120}}
