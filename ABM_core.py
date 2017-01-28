#!/usr/bin/env python
# encoding: utf-8

"""
*****************************************************************************
Author(s):	Maximilian King  (previous authors: Lars Seemann - lseemann@uh.edu)
Email: Maximilian_King@brown.edu
Organization: Marshall Lab, Department of Epidemiology - Brown University

Description:
    Module responsible for ABM simulation events. Operates main loop over simulation run.
    Handles agent pairing, interaction, disease propagation, interventions, deaths, etc.


Copyright (c) 2016, Maximilian King
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************
"""

# Imports
import random
#from copy import deepcopy, copy
import os
import time
#import PyQt4

#import scipy.sparse as spsp
from scipy.stats import binom
from scipy.stats import poisson
from functools import wraps

try:
    from HIVABM_Population import PopulationClass, print_population
except ImportError:
    raise ImportError("Can't import PopulationClass")

try:
    from HIVABM_PartialNetwork import NetworkClass, save_adjlist
except ImportError, e:
    raise ImportError("Can't import NetworkClass! %s" % str(e))

try:
    from partnering import *
except ImportError, e:
    raise ImportError("Can't import assessment_lib! %s" % str(e))

try:
    from analysis_output import * #assessment_lib import *   #OLD FILE
except ImportError, e:
    raise ImportError("Can't import assessment_lib! %s" % str(e))


PROF_DATA = {}
def profile(function):
    @wraps(function)
    def with_profiling(*args, **kwargs):
        start_time = time.time()

        ret = function(*args, **kwargs)

        elapsed_time = time.time() - start_time

        if function.__name__ not in PROF_DATA:
            PROF_DATA[function.__name__] = [0, []]
        PROF_DATA[function.__name__][0] += 1
        PROF_DATA[function.__name__][1].append(elapsed_time)
        #print elapsed_time

        return ret

    return with_profiling

def print_prof_data():
    for fname, data in PROF_DATA.items():
        max_time = max(data[1])
        avg_time = sum(data[1]) / len(data[1])
        print "Function %s called %d times. " % (fname, data[0])
        print '\tExecution time max: %.3f, average: %.3f, total %.3f' % (max_time, avg_time, sum(data[1]))

def clear_prof_data():
    global PROF_DATA
    PROF_DATA = {}


class HIVModel(NetworkClass):
    """
    :Purpose:
        This is the core class used to simulate
        the spread of HIV and drug use in one MSA
        (Metropolitan Statistical Area).

    :Input:
        N : int
            Number of agents. Default: 1000
        tmax: int
            Number of simulation steps (years).

        :py:class:`SocialNetworkClass` : Inherited
        :py:class:`PopulationClass` : Inherited

    :Attributes:
        :py:attr:`tmax` : int
            Number of time steps simulated.

        :py:attr:`CleanSyringeUsers` : list

        :py:attr:`SEPAgents` : dict
            Dictionary of users who participated in a
            syringe exchange program (SEP) {agent:time}.

        :py:attr:`DrugTreatmentAgents` : dict
            Dictionary of users who underwent drug
            treatment {agent:time}.

        :py:attr:`TestedAgents` : list
            List of agents who get tested for HIV every time step.

        :py:attr:`tmp_Agents` : dict
            Changes resulting from parsing through the agents
            and applying the update rules are stored
            in :py:attr:`tmp_agent_dict`.

        All attributes from :py:class:`SocialNetworkClass` \n

        All attributes from :py:class:`PopulationClass`

    :Methods:
        :py:meth:`_update_population` \n
        :py:meth:`_needle_transmission` \n
        :py:meth:`_sex_transmission` \n
        :py:meth:`_drug_transition` \n
        :py:meth:`_update_IDU` \n
        :py:meth:`_update_NIDU_ND` \n
        :py:meth:`_update_AllAgents` \n
        :py:meth:`_VCT` \n
        :py:meth:`_SEP` \n
        :py:meth:`_enter_drug_treatment` \n
        :py:meth:`_initiate_HAART` \n
        :py:meth:'_discontinue_HAART' \n
        :py:meth:`_get_partner` \n
        :py:meth:`_update_AllAgents` \n
        :py:meth:`run` \n
        :py:meth:`store_results` \n
        :py:meth:`get_HIV_prevalence_drugs` \n
        :py:meth:`get_HIV_prevalence` \n
        :py:meth:`plot_results` \n
        All methods from :py:class:`SocialNetworkClass` \n
        All methods from :py:class:`PopulationClass`
    """

    def __init__(self, N, tmax, parameter_dict, rseed, network_type=None):
        """ Initialize HIVModel object """
        if (type(tmax) is not int):
            raise ValueError("Number of time steps must be integer")
        else:
            self.tmax = tmax

        if (type(rseed) is not int):
            raise ValueError("Random seed must be integer")
        else:
            self.randomSeed = rseed

        #random.seed(self.randomSeed)

        self.current_dir = os.getcwd()
        print("\n === Begin Initialization Protocol ===\n")
        self.ExistingLinksCollapsedList = list()

        ############################## NEW PARAMS
        self.PARTNERTURNOVER = .2
        self.HIGH_RISK_PROPORTION = .05
        self.HIGH_RISK_PARTNER_ENHANCEMENT = 3
        ##############################

        self.TIC = 0
        self.TOC = 0

        self.NEEDLESCALINGPARAM = 1.0
        self.SEXSCALINGPARAM = 0.2
        self.XMISSIONSCALINGPARAM = 1.0

        """
        Sensitivity Analysis
        -        Vary infectiousness of acute by 50% and 200% of base value
        -        Risk reduction for PWID DX to fixed 0%, 80%
        -        Risk reduction for PWID HAART to 50% and 200% of base value
        -        Testing frequency for PWID to 50% and 200% of base value
        """

        self.sensAcuteScaling = 1.0
        self.sensRR_Dx = 0.5
        self.sensRR_HAART = 1.0
        self.sensTestFreq = 1.0


        # self.x for combination intervention paperf
        # time dependent vector values, e.g., self.NSP_SAT = dict[time][value]
        self.NSP_NoSAT = parameter_dict['NSP_NoSAT']  # dict
        self.NSP_SAT = parameter_dict['NSP_SAT']  # dict

        # scalar values
        self.SAT_NoNSP_IDU = parameter_dict['SAT_NoNSP_IDU']
        self.SAT_NSP = parameter_dict['SAT_NSP']
        self.SAT_NIDU = parameter_dict['SAT_NIDU']
        self.VCT_NoNSP_IDU = parameter_dict['VCT_NoNSP_IDU']
        self.VCT_NoNSP_MSM = parameter_dict['VCT_NoNSP_MSM']
        self.VCT_NoNSP_EE = parameter_dict['VCT_NoNSP_EE']
        self.VCT_NoNSP_NIDU = parameter_dict['VCT_NoNSP_NIDU']
        self.VCT_NSP = parameter_dict['VCT_NSP']
        self.HAART_NoSAT_IDU = parameter_dict['HAART_NoSAT_IDU']
        self.HAART_NoSAT_NIDU = parameter_dict['HAART_NoSAT_NIDU']
        self.HAART_NoSAT_EE = parameter_dict['HAART_NoSAT_EE']
        self.HAART_SAT_IDU = parameter_dict['HAART_SAT_IDU']
        self.HAART_SAT_NIDU = parameter_dict['HAART_SAT_NIDU']
        self.HAARTdis_SAT_IDU = parameter_dict['HAARTdis_SAT_IDU']
        self.HAARTdis_NoSAT_IDU = parameter_dict['HAARTdis_NoSAT_IDU']
        self.HAARTdis_SAT_NIDU = parameter_dict['HAARTdis_SAT_NIDU']
        self.HAARTdis_NoSAT_NIDU = parameter_dict['HAARTdis_NoSAT_NIDU']
        self.HAARTdis_EE = parameter_dict['HAARTdis_EE']
        self.Adhere_SAT_IDU = parameter_dict['Adhere_SAT_IDU']
        self.Adhere_NoSAT_IDU = parameter_dict['Adhere_NoSAT_IDU']
        self.Adhere_SAT_NIDU = parameter_dict['Adhere_SAT_NIDU']
        self.Adhere_NoSAT_NIDU = parameter_dict['Adhere_NoSAT_NIDU']
        self.Adhere_EE = parameter_dict['Adhere_EE']
        self.SAT_disc = parameter_dict['SAT_disc']

        self.TRMP_NCHRONIC_0 = parameter_dict['TRMP_NCHRONIC_0']
        self.TRMP_NCHRONIC_1 = parameter_dict['TRMP_NCHRONIC_1']
        self.TRMP_NCHRONIC_2 = parameter_dict['TRMP_NCHRONIC_2']
        self.TRMP_NCHRONIC_3 = parameter_dict['TRMP_NCHRONIC_3']
        self.TRMP_NCHRONIC_4 = parameter_dict['TRMP_NCHRONIC_4']
        self.TRMP_NCHRONIC_5 = parameter_dict['TRMP_NCHRONIC_5']

        self.TRMP_NAIDS_0 = parameter_dict['TRMP_NAIDS_0']
        self.TRMP_NAIDS_1 = parameter_dict['TRMP_NAIDS_1']
        self.TRMP_NAIDS_2 = parameter_dict['TRMP_NAIDS_2']
        self.TRMP_NAIDS_3 = parameter_dict['TRMP_NAIDS_3']
        self.TRMP_NAIDS_4 = parameter_dict['TRMP_NAIDS_4']
        self.TRMP_NAIDS_5 = parameter_dict['TRMP_NAIDS_5']

        self.TRMP_NACUTE_0 = parameter_dict['TRMP_NACUTE_0']
        self.TRMP_NACUTE_1 = parameter_dict['TRMP_NACUTE_1']
        self.TRMP_NACUTE_2 = parameter_dict['TRMP_NACUTE_2']
        self.TRMP_NACUTE_3 = parameter_dict['TRMP_NACUTE_3']
        self.TRMP_NACUTE_4 = parameter_dict['TRMP_NACUTE_4']
        self.TRMP_NACUTE_5 = parameter_dict['TRMP_NACUTE_5']

        self.TRMP_CHRONIC_SMSM_0 = parameter_dict['TRMP_CHRONIC_SMSM_0']
        self.TRMP_CHRONIC_SMSM_1 = parameter_dict['TRMP_CHRONIC_SMSM_1']
        self.TRMP_CHRONIC_SMSM_2 = parameter_dict['TRMP_CHRONIC_SMSM_2']
        self.TRMP_CHRONIC_SMSM_3 = parameter_dict['TRMP_CHRONIC_SMSM_3']
        self.TRMP_CHRONIC_SMSM_4 = parameter_dict['TRMP_CHRONIC_SMSM_4']
        self.TRMP_CHRONIC_SMSM_5 = parameter_dict['TRMP_CHRONIC_SMSM_5']

        self.TRMP_CHRONIC_S_0 = parameter_dict['TRMP_CHRONIC_S_0']
        self.TRMP_CHRONIC_S_1 = parameter_dict['TRMP_CHRONIC_S_1']
        self.TRMP_CHRONIC_S_2 = parameter_dict['TRMP_CHRONIC_S_2']
        self.TRMP_CHRONIC_S_3 = parameter_dict['TRMP_CHRONIC_S_3']
        self.TRMP_CHRONIC_S_4 = parameter_dict['TRMP_CHRONIC_S_4']
        self.TRMP_CHRONIC_S_5 = parameter_dict['TRMP_CHRONIC_S_5']

        self.TRMP_AIDS_SMSM_0 = parameter_dict['TRMP_AIDS_SMSM_0']
        self.TRMP_AIDS_SMSM_1 = parameter_dict['TRMP_AIDS_SMSM_1']
        self.TRMP_AIDS_SMSM_2 = parameter_dict['TRMP_AIDS_SMSM_2']
        self.TRMP_AIDS_SMSM_3 = parameter_dict['TRMP_AIDS_SMSM_3']
        self.TRMP_AIDS_SMSM_4 = parameter_dict['TRMP_AIDS_SMSM_4']
        self.TRMP_AIDS_SMSM_5 = parameter_dict['TRMP_AIDS_SMSM_5']

        self.TRMP_AIDS_S_0 = parameter_dict['TRMP_AIDS_S_0']
        self.TRMP_AIDS_S_1 = parameter_dict['TRMP_AIDS_S_1']
        self.TRMP_AIDS_S_2 = parameter_dict['TRMP_AIDS_S_2']
        self.TRMP_AIDS_S_3 = parameter_dict['TRMP_AIDS_S_3']
        self.TRMP_AIDS_S_4 = parameter_dict['TRMP_AIDS_S_4']
        self.TRMP_AIDS_S_5 = parameter_dict['TRMP_AIDS_S_5']

        self.TRMP_ACUTE_SMSM_0 = parameter_dict['TRMP_ACUTE_SMSM_0']
        self.TRMP_ACUTE_SMSM_1 = parameter_dict['TRMP_ACUTE_SMSM_1']
        self.TRMP_ACUTE_SMSM_2 = parameter_dict['TRMP_ACUTE_SMSM_2']
        self.TRMP_ACUTE_SMSM_3 = parameter_dict['TRMP_ACUTE_SMSM_3']
        self.TRMP_ACUTE_SMSM_4 = parameter_dict['TRMP_ACUTE_SMSM_4']
        self.TRMP_ACUTE_SMSM_5 = parameter_dict['TRMP_ACUTE_SMSM_5']

        self.TRMP_ACUTE_S_0 = parameter_dict['TRMP_ACUTE_S_0']
        self.TRMP_ACUTE_S_1 = parameter_dict['TRMP_ACUTE_S_1']
        self.TRMP_ACUTE_S_2 = parameter_dict['TRMP_ACUTE_S_2']
        self.TRMP_ACUTE_S_3 = parameter_dict['TRMP_ACUTE_S_3']
        self.TRMP_ACUTE_S_4 = parameter_dict['TRMP_ACUTE_S_4']
        self.TRMP_ACUTE_S_5 = parameter_dict['TRMP_ACUTE_S_5']

        self.D_IDU_NHIV = parameter_dict['D_IDU_NHIV']
        self.D_IDU_THIV = parameter_dict['D_IDU_THIV']
        self.D_IDU_PHIV = parameter_dict['D_IDU_PHIV']
        self.D_IDU_AIDS = parameter_dict['D_IDU_AIDS']
        self.D_NIDU_NHIV = parameter_dict['D_NIDU_NHIV']
        self.D_NIDU_THIV = parameter_dict['D_NIDU_THIV']
        self.D_NIDU_PHIV = parameter_dict['D_NIDU_PHIV']
        self.D_NIDU_AIDS = parameter_dict['D_NIDU_AIDS']
        self.D_ND_NHIV = parameter_dict['D_ND_NHIV']
        self.D_ND_THIV = parameter_dict['D_ND_THIV']
        self.D_ND_PHIV = parameter_dict['D_ND_PHIV']
        self.D_ND_AIDS = parameter_dict['D_ND_AIDS']

        self.USP_MSM_IDU_UHIV_P1 = parameter_dict['USP_MSM_IDU_UHIV_P1']
        self.USP_MSM_NIDU_UHIV_P1 = parameter_dict['USP_MSM_NIDU_UHIV_P1']
        self.USP_MSM_ND_UHIV_P1 = parameter_dict['USP_MSM_ND_UHIV_P1']
        self.USP_HM_IDU_UHIV_P1 = parameter_dict['USP_HM_IDU_UHIV_P1']
        self.USP_HM_NIDU_UHIV_P1 = parameter_dict['USP_HM_NIDU_UHIV_P1']
        self.USP_HM_ND_UHIV_P1 = parameter_dict['USP_HM_ND_UHIV_P1']
        self.USP_FM_IDU_UHIV_P1 = parameter_dict['USP_FM_IDU_UHIV_P1']
        self.USP_FM_NIDU_UHIV_P1 = parameter_dict['USP_FM_NIDU_UHIV_P1']
        self.USP_FM_ND_UHIV_P1 = parameter_dict['USP_FM_ND_UHIV_P1']
        self.USP_WSW_IDU_UHIV_P1 = parameter_dict['USP_WSW_IDU_UHIV_P1']
        self.USP_WSW_NIDU_UHIV_P1 = parameter_dict['USP_WSW_NIDU_UHIV_P1']
        self.USP_WSW_ND_UHIV_P1 = parameter_dict['USP_WSW_ND_UHIV_P1']

        self.USP_MSM_IDU_UHIV_P2 = parameter_dict['USP_MSM_IDU_UHIV_P2']
        self.USP_MSM_NIDU_UHIV_P2 = parameter_dict['USP_MSM_NIDU_UHIV_P2']
        self.USP_MSM_ND_UHIV_P2 = parameter_dict['USP_MSM_ND_UHIV_P2']
        self.USP_HM_IDU_UHIV_P2 = parameter_dict['USP_HM_IDU_UHIV_P2']
        self.USP_HM_NIDU_UHIV_P2 = parameter_dict['USP_HM_NIDU_UHIV_P2']
        self.USP_HM_ND_UHIV_P2 = parameter_dict['USP_HM_ND_UHIV_P2']
        self.USP_FM_IDU_UHIV_P2 = parameter_dict['USP_FM_IDU_UHIV_P2']
        self.USP_FM_NIDU_UHIV_P2 = parameter_dict['USP_FM_NIDU_UHIV_P2']
        self.USP_FM_ND_UHIV_P2 = parameter_dict['USP_FM_ND_UHIV_P2']
        self.USP_WSW_IDU_UHIV_P2 = parameter_dict['USP_WSW_IDU_UHIV_P2']
        self.USP_WSW_NIDU_UHIV_P2 = parameter_dict['USP_WSW_NIDU_UHIV_P2']
        self.USP_WSW_ND_UHIV_P2 = parameter_dict['USP_WSW_ND_UHIV_P2']

        self.USP_MSM_IDU_PHIV_P1 = parameter_dict['USP_MSM_IDU_PHIV_P1']
        self.USP_MSM_NIDU_PHIV_P1 = parameter_dict['USP_MSM_NIDU_PHIV_P1']
        self.USP_MSM_ND_PHIV_P1 = parameter_dict['USP_MSM_ND_PHIV_P1']
        self.USP_HM_IDU_PHIV_P1 = parameter_dict['USP_HM_IDU_PHIV_P1']
        self.USP_HM_NIDU_PHIV_P1 = parameter_dict['USP_HM_NIDU_PHIV_P1']
        self.USP_HM_ND_PHIV_P1 = parameter_dict['USP_HM_ND_PHIV_P1']
        self.USP_FM_IDU_PHIV_P1 = parameter_dict['USP_FM_IDU_PHIV_P1']
        self.USP_FM_NIDU_PHIV_P1 = parameter_dict['USP_FM_NIDU_PHIV_P1']
        self.USP_FM_ND_PHIV_P1 = parameter_dict['USP_FM_ND_PHIV_P1']
        self.USP_WSW_IDU_PHIV_P1 = parameter_dict['USP_WSW_IDU_PHIV_P1']
        self.USP_WSW_NIDU_PHIV_P1 = parameter_dict['USP_WSW_NIDU_PHIV_P1']
        self.USP_WSW_ND_PHIV_P1 = parameter_dict['USP_WSW_ND_PHIV_P1']

        self.USP_MSM_IDU_PHIV_P2 = parameter_dict['USP_MSM_IDU_PHIV_P2']
        self.USP_MSM_NIDU_PHIV_P2 = parameter_dict['USP_MSM_NIDU_PHIV_P2']
        self.USP_MSM_ND_PHIV_P2 = parameter_dict['USP_MSM_ND_PHIV_P2']
        self.USP_HM_IDU_PHIV_P2 = parameter_dict['USP_HM_IDU_PHIV_P2']
        self.USP_HM_NIDU_PHIV_P2 = parameter_dict['USP_HM_NIDU_PHIV_P2']
        self.USP_HM_ND_PHIV_P2 = parameter_dict['USP_HM_ND_PHIV_P2']
        self.USP_FM_IDU_PHIV_P2 = parameter_dict['USP_FM_IDU_PHIV_P2']
        self.USP_FM_NIDU_PHIV_P2 = parameter_dict['USP_FM_NIDU_PHIV_P2']
        self.USP_FM_ND_PHIV_P2 = parameter_dict['USP_FM_ND_PHIV_P2']
        self.USP_WSW_IDU_PHIV_P2 = parameter_dict['USP_WSW_IDU_PHIV_P2']
        self.USP_WSW_NIDU_PHIV_P2 = parameter_dict['USP_WSW_NIDU_PHIV_P2']
        self.USP_WSW_ND_PHIV_P2 = parameter_dict['USP_WSW_ND_PHIV_P2']

        self.PAIDS_IDU_0 = parameter_dict['PAIDS_IDU_0']
        self.PAIDS_IDU_1 = parameter_dict['PAIDS_IDU_1']
        self.PAIDS_IDU_2 = parameter_dict['PAIDS_IDU_2']
        self.PAIDS_IDU_3 = parameter_dict['PAIDS_IDU_3']
        self.PAIDS_IDU_4 = parameter_dict['PAIDS_IDU_4']
        self.PAIDS_IDU_5 = parameter_dict['PAIDS_IDU_5']

        self.PAIDS_EE_0 = parameter_dict['PAIDS_EE_0']
        self.PAIDS_EE_1 = parameter_dict['PAIDS_EE_1']
        self.PAIDS_EE_2 = parameter_dict['PAIDS_EE_2']
        self.PAIDS_EE_3 = parameter_dict['PAIDS_EE_3']
        self.PAIDS_EE_4 = parameter_dict['PAIDS_EE_4']
        self.PAIDS_EE_5 = parameter_dict['PAIDS_EE_5']

        self.MEAN_N_ACTS = parameter_dict['MEAN_N_ACTS']
        self.MEAN_SMSM_ACTS = parameter_dict['MEAN_SMSM_ACTS']
        self.MEAN_S_ACTS = parameter_dict['MEAN_S_ACTS']

        self.ProbDeath = {'IDU': {}, 'NIDU': {}, 'ND': {}}
        self.ProbDeath['IDU'].update({'HIV-': self.D_IDU_NHIV, 'HIV+/HAART': self.D_IDU_THIV,
                                      'HIV+': self.D_IDU_PHIV, 'AIDS': self.D_IDU_AIDS})
        self.ProbDeath['NIDU'].update({'HIV-': self.D_NIDU_NHIV, 'HIV+/HAART': self.D_NIDU_THIV,
                                       'HIV+': self.D_NIDU_PHIV, 'AIDS': self.D_NIDU_AIDS})
        self.ProbDeath['ND'].update({'HIV-': self.D_ND_NHIV, 'HIV+/HAART': self.D_ND_THIV,
                                     'HIV+': self.D_ND_PHIV, 'AIDS': self.D_ND_AIDS})

        print("\tDictionary Read")

        StratW = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}}
        StratW['MSM'] = {'NUMPartn':1.5, 'NUMSexActs':5.0, 'UNSAFESEX':0.43, 'NEEDLESH':0.43, 'HIVTEST':0.055, 'INCAR':0.00014, 'HAARTadh':0.57, 'HAARTdisc':0.08}
        StratW['HM'] = {'NUMPartn':1.2, 'NUMSexActs':5.0, 'UNSAFESEX':0.83, 'NEEDLESH':0.43, 'HIVTEST':0.013, 'INCAR':0.00014, 'HAARTadh':0.56, 'HAARTdisc':0.01}
        StratW['HF'] = {'NUMPartn':1.2, 'NUMSexActs':5.0, 'UNSAFESEX':0.83, 'NEEDLESH':0.43, 'HIVTEST':0.015, 'INCAR':0.00003, 'HAARTadh':0.53, 'HAARTdisc':0.013}
        StratW['PWID'] = {'NUMPartn':1.2, 'NUMSexActs':5.0, 'UNSAFESEX':0.72, 'NEEDLESH':0.43, 'HIVTEST':0.037, 'INCAR':0.0050, 'HAARTadh':0.52, 'HAARTdisc':0.018}

        StratB = {'MSM':{}, 'HM':{}, 'HF':{}, 'PWID':{}}
        StratB['MSM'] = {'NUMPartn':4.9, 'NUMSexActs':6.0, 'UNSAFESEX':0.43, 'NEEDLESH':0.27, 'HIVTEST':0.06, 'INCAR':0.0011, 'HAARTadh':0.41, 'HAARTdisc':0.01}
        StratB['HM'] = {'NUMPartn':3.5, 'NUMSexActs':5.0, 'UNSAFESEX':0.70, 'NEEDLESH':0.27, 'HIVTEST':0.023, 'INCAR':0.0011, 'HAARTadh':0.34, 'HAARTdisc':0.013}
        StratB['HF'] = {'NUMPartn':3.5, 'NUMSexActs':5.0, 'UNSAFESEX':0.70, 'NEEDLESH':0.27, 'HIVTEST':0.033, 'INCAR':0.00016, 'HAARTadh':0.50, 'HAARTdisc':0.017}
        StratB['PWID'] = {'NUMPartn':3.5, 'NUMSexActs':5.0, 'UNSAFESEX':0.69, 'NEEDLESH':0.27, 'HIVTEST':0.043, 'INCAR':0.0060, 'HAARTadh':0.36, 'HAARTdisc':0.023}

        self.ProbTables = {'WHITE':StratW, 'BLACK':StratB}



        SexTrans = {'MSM':{}, 'HM':{}, 'HF':{}}
        SexTrans['MSM'] = {'0':0.005, '1':0.005, '2':0.004, '3':0.002, '4':0.001, '5':0.0001}
        SexTrans['HM'] = {'0':0.001, '1':0.001, '2':0.0008, '3':0.0004, '4':0.0002, '5':0.0001}
        SexTrans['HF'] = {'0':0.001, '1':0.001, '2':0.0008, '3':0.0004, '4':0.0002, '5':0.0001}

        NeedleTrans = {0:0.007, '1':0.007, '2':0.0056, '3':0.0028, '4':0.0014, '5':0.0002}
        #NeedleTrans = {0.007, 0.007, 0.0056, 0.0028, 0.0014, 0.0002}

        self.TransmissionProbabilities = {'SEX':SexTrans, 'NEEDLE':NeedleTrans}
        # Risk network replaced social network
        if network_type:
            print("\tNetwork Class")
            NetworkClass.__init__(self, N=N, m_0=1, network_type=network_type)
        else:
            print("\tPopulation Class")
            PopulationClass.__init__(self, n=N, rSeed = rseed)

        self.AdjMat = 0#spsp.lil_matrix((self.PopulationSize, self.PopulationSize), dtype=np.int8)
        self.AdjMats_by_time = 0#{0: self.AdjMat}

        # keep track of current time step globally for dynnetwork report
        self.TimeStep = 0


        print("\n\tCreating lists")
        # Other lists / dictionaries
        self.MSM_RandomChoice = []  # List of MSM agents who choose a new
        for agent in self.MSM_agents:
            if agent not in self.IDU_agents:  # partner for interaction randomly
                if random.random() < 0.1:  # MSM agents not in this list choose
                    self.MSM_RandomChoice.append(agent)  # only MSM agents for interactions
        self.MSM_RandomChoice2 = []  # List of MSM - IDU agents who choose a new
        for agent in self.MSM_agents:
            if agent in self.IDU_agents:  # partner for interaction randomly
                if random.random() < 0.5:  # MSM agents not in this list choose
                    self.MSM_RandomChoice2.append(agent)  # only MSM agents for interactions
        self.WSW_RandomChoice = []  # List of WSW agents who choose a new
        for agent in self.WSW_agents:  # partner for interaction randomly
            if random.random() < 0.5:  # WSW agents not in this list choose
                self.WSW_RandomChoice.append(agent)  # only WSW agents for interactions
        self.SEPAgents = dict()  # dict of users who used SEP (agent:time)
        self.SEPAgents_past = dict()  # dict of users who used SEP in last time step (agent:time)
        self.N_TreatmentSpots = 0  # initiate number of occupied treatment spots

        self.DrugTreatmentAgents_current = dict()  # dictionary of users who are currently undergoing
        for agent in self.IDU_agents or self.NIDU_agents:  # drug treatent (agent:time)
            if random.random() < 0.09:  # Initialize
                self.DrugTreatmentAgents_current.update({agent: 1})
        self.DrugTreatmentAgents_past = dict()  # dictionary of users who underwent drug treat-
        # ment in the past(agent:time)
        self.VCTAgents = dict()  # list of agents who get tested for HIV ((agent:time)
        # dict of agents HIV treatment (HAART) adherence (agent:adherence)
        self.AdherenceAgents = dict()
        for agent in self.Agents:
            self.AdherenceAgents.update({agent: 0})
        self.Viral_load = dict()
        for agent in self.Agents:
            self.Viral_load.update({agent: 0})
        # dict of agents who have HIV in the past time step
        self.HIV_key_transitiontime = {1: copy(self.HIV_agents)}
        for t in range(2, tmax + 1):
            self.HIV_key_transitiontime.update({t: []})

        self.high_risk_agents = set()
        for agent in self.Agents:
            if np.random.uniform() < self.HIGH_RISK_PROPORTION:
                self.high_risk_agents.add(agent)

        self.NewInfections = Agent_set(3, "NewInfections")
        # Tmp agent lists to keep track of updates
        # self.tmp_Agents = deepcopy(self.Agents)  # dictionary that keeps track of the update
        #
        # self.tmp_IDU_agents = copy(self.IDU_agents)
        # self.tmp_NIDU_agents = copy(self.NIDU_agents)
        # self.tmp_ND_agents = copy(self.ND_agents)
        #
        # self.tmp_HM_agents = copy(self.HM_agents)
        # self.tmp_HF_agents = copy(self.HF_agents)
        # self.tmp_MSM_agents = copy(self.MSM_agents)
        # self.tmp_WSW_agents = copy(self.WSW_agents)
        #
        # self.tmp_HIV_agents = copy(self.HIV_agents)
        # self.tmp_AIDS_agents = copy(self.AIDS_agents)
        # self.tmp_HAART_agents = copy(self.HAART_agents)

        # Keep track of HIV results
        self.HIVinIDU = [0] * self.tmax
        self.HIVinNIDU = [0] * self.tmax
        self.HIVinND = [0] * self.tmax
        self.HIVinMSM = [0] * self.tmax
        self.HIVinHM = [0] * self.tmax
        self.HIVinHF = [0] * self.tmax
        self.HIVinWSW = [0] * self.tmax
        self.HIVinMIDU = [0] * self.tmax
        self.HIVinMNIDU = [0] * self.tmax
        self.HIVinMND = [0] * self.tmax
        self.HIVinIDUnmsm = [0] * self.tmax
        self.HIVinNIDUnmsm = [0] * self.tmax
        self.HIVinNDnmsm = [0] * self.tmax

        # Assess the distribution of number of interactions per timestep for each agent type
        self.ND_NumPartners = {'ND': [], 'NIDU': [], 'IDU': [], 'MSM': []}  # final counts
        self.NIDU_NumPartners = {'ND': [], 'NIDU': [], 'IDU': [], 'MSM': []}  # final counts
        self.IDU_NumPartners = {'ND': [], 'NIDU': [], 'IDU': [], 'MSM': []}  # final counts
        self.MSM_NumPartners = {'ND': [], 'NIDU': [], 'IDU': [], 'MSM': []}  # final counts

        self.tmp_ND_NumPartners_Count = {}
        self.tmp_NIDU_NumPartners_Count = {}
        self.tmp_IDU_NumPartners_Count = {}
        self.tmp_MSM_NumPartners_Count = {}
        self.tmp_WSW_NumPartners_Count = {}

        self.Acute_agents = []
        self.Transmit_from_agents = []
        self.Transmit_to_agents = []
        self.Transmission_tracker = {'SEX_MSM': {1: 0}, 'SEX_NMSM': {1: 0}, 'NEEDLE': {1: 0}}

        self.ResultDict = initiate_ResultDict()

        print("\tReseting death count")
        self._reset_death_count()  # Number of death
        print("\n === Initialization Protocol Finished === \n")

    def run(self, save_adjlist_flag=1, dir_prefix='Results'):
        """
        Core of the model:
            1. Prints networkReport for first agents.
            2. Makes agents become HIV (used for current key_time tracking for acute)
            3. Loops over all time steps
                a. _update AllAgents()
                b. _reset_death_counts()
                c. _ self._die_and_replace()
                d. self._update_population()
                e. self._reset_partner_count()
        """

        #print "RANDOM CALL %d" %random.randint(0,100)
        #random.seed(self.randomSeed)

        #print "RANDOM CALL %d" %random.randint(0,100)
        #@profile
        def getStats(t):
            self.filler = 0
            print_stats(t, self.totalAgentClass, self.HIV_agents_class, self.IncarceratedClass, self.PrEP_agents_class, self.NewInfections, self.num_Deaths)

            """self.ResultDict = assess_before_update(t,
                                                   self.ResultDict,
                                                   self.Agents,
                                                   self.HIV_agents,
                                                   self.filler,
                                                   self.AIDS_agents,
                                                   self.filler,
                                                   self.SEPAgents,
                                                   self.HAART_agents,
                                                   self.filler,
                                                   self.AdherenceAgents,
                                                   self.num_Deaths,
                                                   self.AdjMats_by_time,
                                                   self.Acute_agents,
                                                   self.Transmit_from_agents,
                                                   self.Transmit_to_agents,
                                                   self.Transmission_tracker,
                                                   self.high_risk_agents,
                                                   self.Incarcerated,
                                                   self.HIVidentified_agents)
            """

        #random.seed(self.randomSeed)

        print("\n === Begin Simulation Run ===\n")
        print("\t Writing Agents to dynNet Report")

        # write agents to dynnetworkReport
        """
        dynnetworkReport = open('Results/dynnetworkReport.txt', 'a')
        for agent in self.Agents:
            sextype = self.Agents[agent]['Sex Type']
            drugtype = self.Agents[agent]['Drug Type']
            HIV = self.Agents[agent]['HIV']
            reportLine = '\t'.join(['0', 'NEWAGENT', repr(agent), sextype, drugtype, repr(HIV)])
            dynnetworkReport.write('\n' + reportLine)
            #open('dynnetworkReport.txt', 'a').write('\n' + reportLine)

        dynnetworkReport.close()
        # initiate HIV status
        for agent in self.Agents:
            if self.Agents[agent]['HIV'] == 1: self._become_HIV(agent, 0)
        #HIVcountReport = open('Results/HIV_Report.txt', 'w')



        startReport = open('Results/startingAgentReport.txt', 'w')
        startReport.write('Agent\tRace\tSex\tDrug\tHIV\tAIDS\tHAART\tIncar\tIncar_t\n')
        for agent in self.Agents:
            racetype = self.Agents[agent]['Race']
            sextype = self.Agents[agent]['Sex Type']
            drugtype = self.Agents[agent]['Drug Type']

            HIV = self.Agents[agent]['HIV']
            AIDS = self.Agents[agent]['AIDS']
            HAART = self.Agents[agent]['HAARTa']
            incar_t = self.Agents[agent]['incar_t']
            if incar_t > 0:
                incar = 1
            else:
                incar = 0

            if HAART == 1: ### FIXES ADHERENCE FOR INITIALIZED AGENTS
                self._initiate_HAART(agent, 0)
            #reportLine = agent, racetype,sextype,drugtype,HIV, AIDS, HAART,incar
            startReport.write('%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n'%(agent, racetype,sextype,drugtype,HIV, AIDS, HAART,incar, incar_t))
            #open('Results/dynnetworkReport.txt', 'a').write('\n' + reportLine)

        startReport.close()
        """

        print "TOTAL AGENTS:", self.totalAgentClass.num_members()

        #self.totalAgentClass.print_agents_to_file(0,"wb")
        #self.NewInfections.print_agents_to_file(0,"wb","NewInfections.txt")
        #self.Relationships.print_agent_relationships_to_file(0,"wb")
        print("\t===! Start Main Loop !===\n")
        for t in range(1, self.tmax + 1):
            print '\n\t\t\t\t\t\t\t\t\t\t\t\t\t.: TIME', t
            print "RANDOM CALL %d" %random.randint(0,100)
            hivCount = len(self.HIV_agents)

            #todo: GET THIS TO THE NEW HIV COUNT
            #print "\t\tStarting HIV count:", hivCount
            print "\t\tSTARTING HIV count:%d\tIncarcerated:%d\tHR+:%d" % (
            self.totalAgentClass._subset["HIV"].num_members(), self.IncarceratedClass.num_members(), self.PrEP_agents_class.num_members())
            #self.totalAgentClass.print_agents()
            self.TimeStep = t
            self.HIV_key_transitiontime[t] = []

            self._update_AllAgents(t)


            print "Copying ADJMats"
            #self.AdjMats_by_time.update({t: deepcopy(self.AdjMat)})

            print "Transmission Tracker"
            for k in self.Transmission_tracker: self.Transmission_tracker[k].update({t + 1: 0})

            print "Results Dictionary update"
            getStats(t)

            print "Reseting death count"
            self._reset_death_count()

            print("\t\tdie and replace")
            self._die_and_replace()

            #print("\t\tupdate population")
            self._update_population()

            #print("\t\tassess after population update")
            self._reset_partner_count()

            print "\t\tENDING HIV count:%d\tIncarcerated:%d\tHR+:%d"%(self.totalAgentClass._subset["HIV"].num_members(), self.IncarceratedClass.num_members(),self.PrEP_agents_class.num_members())
            #self.Relationships.print_agent_relationshps()
            #print "RANDOM CALL %d\n" %random.randint(0,100)

            #self.totalAgentClass.print_agents_to_file(t)
            #self.Relationships.print_agent_relationships_to_file(t)
            self.totalAgentClass.print_subsets()
            #self.NewInfections.print_agents_to_file(t,filename="NewInfections.txt")
            #self.NewInfections.print_agents()
            self.NewInfections.clear_set()
            self.num_Deaths
            #self.NewInfections.print_agents()
            # print "MSM AGENTS", self.MSM_agentsClass.num_members()
            # self.MSM_agentsClass.print_agents()
            # print "HM AGENTS", self.HM_agentsClass.num_members()
            # self.HM_agentsClass.print_agents()
            # print "HF AGENTS:", self.HF_agentsClass.num_members()
            # self.HF_agentsClass.print_agents()
        #self.totalAgentClass.print_agents_to_file(t)
        endReport = open('Results/endingAgentReport.txt', 'w')
        endReport.write('Agent\tRace\tSex\tDrug\tHIV\tAIDS\tHAART\tIncar\n')
        for agent in self.Agents:
            racetype = self.Agents[agent]['Race']
            sextype = self.Agents[agent]['Sex Type']
            drugtype = self.Agents[agent]['Drug Type']

            HIV = self.Agents[agent]['HIV']
            AIDS = self.Agents[agent]['AIDS']
            HAART = self.Agents[agent]['HAARTa']
            incar_t = self.Agents[agent]['incar_t']
            if incar_t > 0:
                incar = 1
            else:
                incar = 0
            #reportLine = agent, racetype,sextype,drugtype,HIV, AIDS, HAART,incar
            endReport.write('%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n'%(agent, racetype,sextype,drugtype,HIV, AIDS, HAART,incar,incar_t))
            #open('Results/dynnetworkReport.txt', 'a').write('\n' + reportLine)

        endReport.close()
        print_prof_data()

        #self.totalAgentClass.print_agents()



    #@profile
    def _update_AllAgents(self, time):
        """
        :Purpose:
            Update IDU agents:
            For each agent:
                1 - determine agent type
                2 - get partners
                3 - agent interacts with partners
                4 - drug transition
                5 - VCT (Voluntsry Counseling and Testing)
                6 - if IDU: SEP, treatment
                7 - if HIV: HAART, AIDS
                8 - drug cessation

        :Input:
            agent, time

        :Output:
            none
®
        """
        num_HIV = len(self.HIV_agents)
        print("\t\t= Begin Agents Partnering =")
        if time == 1:
            update_partner_assignments(self, 10.0)
            #update_partner_assignments(self, self.PARTNERTURNOVER)
        else:
            update_partner_assignments(self, self.PARTNERTURNOVER)
        print("\t\t= Updated Partners =")
        self.Acute_agents = []
        self.Transmit_from_agents = []
        self.Transmit_to_agents = []
        # self.totalAgentClass.print_agents()
        # write partnerships to dynnetworkReport
        """
        dynnetworkReport = open('Results/dynnetworkReport.txt', 'a')
        for agent in self.Agents:
            for partner in list(self.AdjMat.rows[agent]):
                reportLine = '\t'.join([repr(time), 'PARTNERSHIP', repr(agent), repr(partner)])
                #open('Results/dynnetworkReport.txt', 'a').write('\n' + reportLine)
                dynnetworkReport.write('\n' + reportLine)
        dynnetworkReport.close()
        """
        #print("\t\t= Printed DynNetwork Report =")
        print("\n\t\t= Begin Agents Operations =")
        #todo: Test this new structure
        print("\n\t\t= Relationship Iterations =")
        for rel in self.Relationships.iter_agents():
            self._agents_interact(rel._ID1, rel._ID2, time, rel)
            #print self.Relationships.num_members()
            if rel.progress(1):
                self.Relationships.remove_agent(rel)
            #print self.Relationships.num_members()

        for tmpA in self.HighriskClass.iter_agents():
            if tmpA._highrisk_time > 0:
                tmpA._highrisk_time -= 1
            else:
                self.HighriskClass.remove_agent(tmpA)
                tmpA._highrisk_bool = False
                if tmpA._SO == "HM":
                    tmpA._mean_num_partners = 3
                elif tmpA._SO == "HF":
                    tmpA._mean_num_partners = 2

        print("\n\t\t= Agents Iterations (Incar/test/AIDS/HAART/PrEP =")
        #random.shuffle(self.Agents)    ###########ADDED SHUFFLE AGENT LIST
        for agent in self.totalAgentClass.iter_agents():#self.Agents: #NEW METHOD
            #print("\nAgent:",agent)
            # agent_dict = self.Agents[agent]
            agent_drug_type = agent._DU# agent_dict['Drug Type']
            agent_sex_type = agent._SO#agent_dict['Sex Type']
            agent_HIV_status = agent._HIV_bool#agent_dict['HIV']

            agent_incarcerated = agent._incar_bool#agent_dict['incar_t']

            # Agent - Partner interaction
            partners = agent._partners#list(self.AdjMat.rows[agent])
            #print agent_HIV_status
            #print("\tBegin partner interactions:")
            if partners and agent_HIV_status and agent_incarcerated is False:
                #print "Agent %d has partners:" % agent._ID, partners
                #intersection = list(set(need_new_partners[:100]).intersection(set(x)))
                #HIV_partners = list(set(self.AdjMat.rows[agent])).intersection(set(self.tmp_HIV_agents))
                #print HIV_partners
                """for num_Interaction, partner in enumerate(partners):
                    #print "\tAgent %d is now interacting with partner %d" % (agent, partner)
                    self._agents_interact(agent, partner, time, num_Interaction)
                    #self._drug_transition(agent, partner)"""
                for rel in agent._relationships:
                    #print "\tAgent %d is now interacting with partner %d" % (rel._ID1.get_ID(), rel._ID2.get_ID())
                    #self._agents_interact(rel._ID1, rel._ID2, time, 10)
                    #rel.progress(1)
                    pass
            #print("\tEnded partner interactions\n")
            #self._VCT(agent, time)  # VCT: Counseling and Testing
            #print "\tIncarcerate Agents"

            if agent._gender == "Q":
                self._incarcerate(agent, time) ##########*****************************************

            #if agent_drug_type == 'IDU':
                #self._SEP(agent, time)  # SEP: Syringe Exchange program
            """########################################PUT THESE BACK IN
            if agent_drug_type in ['NIDU', 'IDU']:
                #print("\tDrug Cessation")
                self._drug_cessation(agent, agent_drug_type)
                #print("\tEnter/Exit Drug Treatment")
                self._enter_and_exit_drug_treatment(agent, time)
            """

            if agent_HIV_status:
                #print "Viral Load Reset"
                #self._become_HIV(agent, time)  # reset viral load ############# TURNED OFF NO NEED TO RESET VIRAL LOAD
                #P_HAART = random.random()
                #if P_HAART < 0.005:
                self._drugTest(agent, time)
                #self._initiate_HAART(agent, time)
                self._progress_to_AIDS(agent, agent_drug_type)
                #hiv_t = agent_dict['HIV_t']
                #hiv_t += 1
                agent._HIV_time += 1
                # self.tmp_Agents[agent].update({'HIV_t': hiv_t})
            else:
                if not agent._PrEP_bool:
                    #pass
                    self._initiate_PrEP(agent, time)
        print("\t\t= End Agents Operations =")

        """
        # Check consistency
        #print("\n\t\t= Check Consistency =")
        if len(self.Agents) != self.PopulationSize:
            raise ValueError("Check population size! %d - %d" % (
                len(self.Agents), self.PopulationSize))
        if len(self.HIV_agents) != num_HIV:
            raise ValueError("self.HIV_agents modified in _update_AllAgents!")

        """
        print("\t\t !!!! ALL AGENTS UPDATED !!!\n")


    def _agents_interact(self, agent, partner, time, rel):
        """
        :Purpose:
            Let IDU agent interact with a partner.
            Update IDU agents:
                1 - determine transition type
                2 - Injection rules
                3 - Sex rules
                4 - HIV transmission
                5 - SEP

        :Input:
            agent : int

            partner : int

            time : int

        Output:
            none

        """
        #print agent
        partner_drug_type = partner._DU#self.get_agent_characteristic(partner, 'Drug Type')
        agent_drug_type = agent._DU#self.get_agent_characteristic(agent, 'Drug Type')
        partner_sex_type = partner._SO#self.get_agent_characteristic(partner, 'Sex Type')
        agent_sex_type = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        partner_HIV_status = partner._HIV_bool#self.get_agent_characteristic(partner, 'HIV')
        agent_HIV_status = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        agent_incar = agent._incar_bool
        partner_incar = partner._incar_bool
        num_interactions = random.randrange(1, 50, 1)
        if agent_HIV_status == 1 and partner_HIV_status == 1:
            #print "\t\t!!! Both agents HIV+, moving on A:%d P:%d"%(agent, partner)
            return
        elif agent_HIV_status == 0 and partner_HIV_status == 0:
            #print "\t\t!!! Neither agents HIV+, moving on A:%d P:%d"%(agent, partner)
            return
        elif agent_incar or partner_incar:
            return
        else:
            if partner_drug_type == 'IDU' and agent_drug_type == 'IDU':
                # Injection is possible

                if self._sex_possible(agent_sex_type, partner_sex_type):
                    # Sex is possible
                    rv = random.random()
                    if rv < 0.6: #Needle only (60%)
                        #print "Needle inc (IDUs)"
                        self._needle_transmission(agent, partner, time)
                    elif rv < 0.6 + 0.2: #Sex only (20%)
                        #print "Sex inc (IDUs)"
                        self._sex_transmission(agent, partner, time)  # , num_interactions)
                    else: #Both sex and needle (20%)
                        #print "Needle and sex inc (IDUs)"
                        self._needle_transmission(agent, partner, time)
                        self._sex_transmission(agent, partner, time)  # , num_interactions)
                else:
                    # Sex not possible, needle only
                    #print "Needle inc (IDUs)"
                    self._needle_transmission(agent, partner, time)

            elif (partner_drug_type in ['NIDU', 'ND'] or agent_drug_type in ['NIDU', 'ND']):
                #print "Sex inc (ND/NIDU)"
                if self._sex_possible(agent_sex_type, partner_sex_type):
                    self._sex_transmission(agent, partner, time, rel)  # ,num_interactions)
                else:
                    return
                    #print "!!!!!!!!!!!!SEX NOT POSSIBLE A:%s \tP:%s"%(agent_sex_type, partner_sex_type)
            else:
                raise ValueError("Agents must be either IDU, NIDU, or ND")



    def _drug_transition(self, agent, partner):
        """
        :Purpose:
            Simulate transition of drug behavior. The following scenarios are
            possible:
            + ND agent might become NIDU when meeting NIDU
            + NIDU might become IDU when meeting IDU
            The function is only applied for NIDU and ND users.

        :Input:
            agents : int
            partner : int

        :Output: -
        """

        partner_drug_type = self.get_agent_characteristic(partner, 'Drug Type')
        agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
        Flag_Partner_IDU_NIDU_Transition = 0
        Flag_Agent_IDU_NIDU_Transition = 0

        # NIDU -> IDU
        if agent_drug_type == 'NIDU' and partner_drug_type == 'IDU':
            if random.random() < 0.00875 / 12:
                #print "Agent %d just became IDU" % (agent, )
                self.tmp_Agents[agent].update({'Drug Type': 'IDU'})  # agent becomes IDU
            # Sex type lists
            if agent in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.remove(agent)
            if agent in self.tmp_ND_agents:  # agent might have transitioned into ND before
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(agent)

        elif partner_drug_type == 'NIDU' and agent_drug_type == 'IDU':
            if random.random() < 0.0175 / 12:
                self.tmp_Agents[partner].update({'Drug Type': 'IDU'})  # partner becomes IDU
            # Sex type lists
            if partner in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.remove(partner)
            if partner in self.tmp_ND_agents:  # agent might have transitioned into ND before
                self.tmp_ND_agents.remove(partner)
            if partner not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(partner)

        ## ND -> IDU
        elif agent_drug_type == 'ND' and partner_drug_type == 'IDU':
            if random.random() < 0.001:
                self.tmp_Agents[agent].update({'Drug Type': 'IDU'})  # agent becomes NIDU
            # Sex type lists
            if agent in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(agent)

        # ND -> IDU
        elif partner_drug_type == 'ND' and agent_drug_type == 'IDU':
            if random.random() < 0.001:
                self.tmp_Agents[agent].update({'Drug Type': 'IDU'})  # agent becomes NIDU
            # Sex type lists
            if agent in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_IDU_agents:
                self.tmp_IDU_agents.append(agent)

        # ND -> NIDU
        elif agent_drug_type == 'ND' and partner_drug_type == 'NIDU':
            if random.random() < 0.005:
                self.tmp_Agents[agent].update({'Drug Type': 'NIDU'})  # agent becomes NIDU
            # Sex type lists
            if agent in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(agent)
            if agent not in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.append(agent)
        # ND -> NIDU
        elif partner_drug_type == 'ND' and agent_drug_type == 'NIDU':
            if random.random() < 0.005:
                self.tmp_Agents[partner].update({'Drug Type': 'NIDU'})  # partner becomes NIDU
            # Sex type lists
            if partner in self.tmp_ND_agents:  # agent might have transitioned already
                self.tmp_ND_agents.remove(partner)
            if partner not in self.tmp_NIDU_agents:
                self.tmp_NIDU_agents.append(partner)
                # if partner in self.tmp_IDU_agents:     # agent might have previously transitioned into IDU
            #    self.tmp_IDU_agents.remove(agent)

        # NIDU -> ND from agent's perspective
        elif agent_drug_type == 'NIDU' and partner_drug_type == 'ND':
            if random.random() < 0.001:
                self.tmp_Agents[agent].update({'Drug Type': 'ND'})
            # Sex type lists
            if agent in self.tmp_NIDU_agents:  # agent might have transitioned already
                self.tmp_NIDU_agents.remove(agent)
            if agent in self.tmp_IDU_agents:  # agent might have previously transitioned into IDU
                self.tmp_IDU_agents.remove(agent)
            if agent not in self.tmp_ND_agents:
                self.tmp_ND_agents.append(agent)
            if agent in self.DrugTreatmentAgents_current:
                self._exit_drug_treatment(agent)

        # NIDU -> ND from partner's perspective
        elif partner_drug_type == 'NIDU' and agent_drug_type == 'ND':
            if random.random() < 0.001:
                self.tmp_Agents[partner].update({'Drug Type': 'ND'})
            # Sex type lists
            if partner in self.tmp_NIDU_agents:  # partner might have transitioned already
                self.tmp_NIDU_agents.remove(partner)
            if partner in self.tmp_IDU_agents:  # partner might have previously transitioned into IDU
                self.tmp_IDU_agents.remove(partner)
            if partner not in self.tmp_ND_agents:
                self.tmp_ND_agents.append(partner)
            if partner in self.DrugTreatmentAgents_current:
                self._exit_drug_treatment(partner)
        else:
            pass  # transition not possible


    def get_acute_status(self, agent, time):
        """
        :Purpose:
            Simulate random transmission of HIV between two IDU agents
            through needle.\n
            Needed in _update_IDUand
        :Input:
            agents : int
            partner : int
        time : int
        :Output: -
        """
        acuteTimePeriod = 3
        hiv_t = agent._HIV_time#self.Agents[agent]['HIV_t']



        #if time > 0:
            #print self.HIV_key_transitiontime
        if hiv_t <= acuteTimePeriod and hiv_t > 0:
            #print "Agent %d has been sick for %d timesteps"%(agent, hiv_t)
            return True
        else:
            return False


    def get_transmission_probability(self, agent, interaction):
        """ Decriptor
            :Purpose:
            Determines the probability of a transmission event based on type. Determines if act is needle/sexual,

            :Input:
                N : int
                Number of agents. Default: 1000
                tmax: int
                Number of simulation steps (years).

                :py:class:`SocialNetworkClass` : Inherited
                :py:class:`PopulationClass` : Inherited

            :Attributes:
                :py:attr:`tmax` : int
                    Number of time steps simulated.
                """

        sex_type = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        race_type = agent._race#self.get_agent_characteristic(agent, 'Race')
        #time = self.TimeStep
        tested = agent._tested#self.get_agent_characteristic(agent, 'Tested')
        onHAART = agent._HAART_bool
        # viral_load = self.Viral_load[agent]
        # v = 10 ** viral_load
        # p = (1. - (1. - (0.317 * (v ** 1.02)) / (v ** 1.02 + 13938 ** 1.02)) ** (1. / 83.17544)) * (0.014 / 0.003)
        #print "VIRAL LOAD: %.2lf, p: %.5lf"%(viral_load, p)
        #isAcute = False
        agentAdherence = agent._HAART_adh#str(self.AdherenceAgents[agent])
        "Logic for if needle or sex type interaction"
        if interaction == 'NEEDLE':
            #p = 0.007
            p = self.TransmissionProbabilities['NEEDLE'][agentAdherence]
            #print p
            #p = p   ################ OLD style using viral load

        elif interaction == 'SEX':
            # agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
            #print sex_type, agentAdherence
            p = self.TransmissionProbabilities['SEX'][sex_type][str(agentAdherence)]
            #p=1#p = 0.0001
            #print p
            """
            if self.get_agent_characteristic(agent, 'Sex Type') == 'MSM':
                #p = 0.005
                p = p * .75  # 3/4
            else:
                #p = 0.001
                p = p / (0.014 / 0.003)  # ~1/4
            """

        isAcute = self.get_acute_status(agent, 0)

        if(isAcute):
            pass
            #p = p * 10 * self.sensAcuteScaling
            #HIV_agent = self.get_agent_characteristic(agent, 'HIV')
            #print "Acute agent probability (%d) @ %.5lf %d" % (agent, p, HIV_agent)

        if tested:
            p = p * 0.5#- (p * self.sensRR_Dx)
            #print "POSITIVE TEST REDUCTION"

        if onHAART:#self.AdherenceAgents[agent] > 0:
            p = p * self.sensRR_HAART


        p = p * self.XMISSIONSCALINGPARAM
        #p = 0
        return p

    #@profile
    def _needle_transmission(self, agent, partner, time):
        """
        :Purpose:
            Simulate random transmission of HIV between two IDU agents
            through needle.\n
            Needed in _update_IDUand
        :Input:
            agents : int
            partner : int
        time : int
        :Output: -
        """

        #Param to scale number of partners
        #NEEDLESCALINGPARAM = 0.7

        # both must be IDU
        partner_drug_type = partner._DU#self.get_agent_characteristic(partner, 'Drug Type')
        agent_drug_type = agent._DU#self.get_agent_characteristic(agent, 'Drug Type')
        agent_race = agent._race#self.get_agent_characteristic(agent, 'Race')
        agent_sex_type = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        Race_Agent = agent._race#self.get_agent_characteristic(agent, 'Race')
        Type_agent = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')


        if not (partner_drug_type == 'IDU' and agent_drug_type == 'IDU'):
            raise ValueError("To share a needle both agents must be IDU!%s %s" %
                             (str(agent_drug_type), str(partner_drug_type)))
        NumberP = len(self.ExistingLinksCollapsedList)
        # Do they share a needle?
        # OLD: if (agent in self.SEPAgents or partner in self.SEPAgents):
        #SEPstat = self._SEP(agent, time)
        SEPstat = False

        """
        tmpRand = random.random()
        if (tmpRand > 0.43 and agent_race == 'WHITE') or (tmpRand > 0.27 and agent_race == 'BLACK'):################ FIX FOR RACIAL DISPARITY
            #print "SEP SAVE %s %d rolled %.3lf"%(agent_race, agent, tmpRand)
            SEPstat = True

        """

        isAcute = self.get_acute_status(agent, time)

        if SEPstat:
            #print "SEP SAVE"
            pass  # no needle sharing
        else:  # they do share a needle
            # HIV+ ?
            HIV_agent = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
            HIV_partner = partner._HIV_bool#self.get_agent_characteristic(partner, 'HIV')
            MEAN_N_ACTS = self.ProbTables[Race_Agent][Type_agent]['NUMSexActs'] * self.NEEDLESCALINGPARAM
            share_acts = poisson.rvs(MEAN_N_ACTS, size=1)
            #share_acts = int(random.uniform(1,30))
            if share_acts < 1:
                share_acts = 1

            p_UnsafeNeedleShare = self.ProbTables[agent_race][agent_sex_type]['NEEDLESH']
            #MSexActs = self.ProbTables[Race_Agent][Type_agent]['NUMSexActs']
            for n in range(share_acts):
                if random.random() > p_UnsafeNeedleShare:
                    share_acts -= 1

            if HIV_agent == 1 and HIV_partner == 0 and share_acts >= 1.0:
                p = self.get_transmission_probability(agent, 'NEEDLE')
                #print p
                p_transmission = binom.pmf(1.0, share_acts, p)

                p_total_transmission = 0
                if share_acts == 1:
                    p_total_transmission = p
                else:
                    for k in range(1, share_acts+1):
                        temp = binom.pmf(k, share_acts, p)
                        #print temp
                        p_total_transmission += temp


                #print "\t\t\tHIV+ NED Act A:%d on P:%d\t Must be less than %.10lf\t(k:1 n:%.2lf, p:%.5lf) **OLD pTrans:%.5lf\tAcute:%s" % (agent, partner, p_total_transmission, share_acts, p, p_transmission, isAcute)
                if random.random() < p_total_transmission:
                    # if agent HIV+ partner becomes HIV+
                    #transmit_HIV(agent, partners, t, )
                    self._become_HIV(partner, time)
                    self.Transmission_tracker['NEEDLE'][time] += 1
                    self.Transmit_from_agents += [agent]
                    self.Transmit_to_agents += [partner]
                    #if agent in ((self.HIV_key_transitiontime[time - 1] if time > 1 else [])
                    #                 + (self.HIV_key_transitiontime[time - 2] if time > 2 else [])
                    #                 + (self.HIV_key_transitiontime[time - 3] if time > 3 else [])):
                    if(isAcute):
                        self.Acute_agents += [agent]
                        print "\t\t(ACUTE)\tNE_HIV from agent %d to partner %d \t@ p=%.5lf transmissionp=%.5lf n:%d" %(agent._ID, partner._ID, p, p_total_transmission, share_acts)
                    else:
                        print "\t\t\t\tNE_HIV from agent %d to partner %d \t@ p=%.5lf transmissionp=%.5lf n:%d" %(agent._ID, partner._ID, p, p_total_transmission, share_acts)

            """elif HIV_partner == 1 and HIV_agent == 0:
                p = self.get_transmission_probability(partner, 'NEEDLE')
                p_transmission = binom.pmf(1.0, share_acts, p)
                print "\t\t\tHIV+ NED Act P:%d on A:%d\t Must be less than %.10lf\t(k:1 n::%.2lf, p:%.5lf)" % (partner, agent, p_transmission, share_acts, p)
                if random.random() < p_transmission:
                    # if partner HIV+ agent becomes HIV+
                    #print "\t\t\t\tNE_HIV from partner %d to agent %d @ p=%.5lf transmissionp=%.5lf" %(partner, agent, p, p_transmission)
                    self._become_HIV(agent, time)
                    self.Transmission_tracker['NEEDLE'][time] += 1
                    self.Transmit_from_agents += [partner]
                    #if partner in ((self.HIV_key_transitiontime[time - 1] if time > 1 else [])
                    #                   + (self.HIV_key_transitiontime[time - 2] if time > 2 else [])
                    #                   + (self.HIV_key_transitiontime[time - 3] if time > 3 else [])):
                        #self.Acute_agents += [partner]
                        #print '\t\t\t\tACUTE', agent
                    if(self.get_acute_status(partner, time)):
                        self.Acute_agents += [partner]
                        print "\t\t\t\tNE_HIV (ACUTE) from partner %d to agent %d \t@ p=%.5lf transmissionp=%.5lf" %(partner, agent, p, p_transmission)
                    else:
                        print "\t\t\t\tNE_HIV from partner %d to agent %d \t@ p=%.5lf transmissionp=%.5lf" %(partner, agent, p, p_transmission)"""

    #@profile
    def _sex_transmission(self, agent, partner, time, rel):

        """
        :Purpose:
            Simulate random transmission of HIV between two agents through Sex.
            Needed for all users. Sex is not possible in case the agent and
            assigned partner have incompatible Sex behavior.

        :Input:
            agents : int
            partner : int
            time : int
        number_of_interaction : int

        :Output:
            none
        """
        #SEXSCALINGPARAM = 0.2

        # Double check: Sex possible?
        Type_agent = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        Type_partner = partner._SO#self.get_agent_characteristic(partner, 'Sex Type')
        if not self._sex_possible(Type_agent, Type_partner):
            raise ValueError("Sex must be possible! %s %s" % (
                str(Type_agent), str(Type_partner)))
        #NumberP = len(self.ExistingLinksCollapsedList)

        # HIV status of agent and partner
        # Everything from here is only run if one of them is HIV+
        HIVstatus_Agent = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        HIVstatus_Partner = partner._HIV_bool#self.get_agent_characteristic(partner, 'HIV')
        AIDSstatus_Agent = agent._AIDS_bool#self.get_agent_characteristic(agent, 'AIDS')
        AIDSstatus_Partner = partner._AIDS_bool#self.get_agent_characteristic(partner, 'AIDS')
        Race_Agent = agent._race#self.get_agent_characteristic(agent, 'Race')
        isAcute = self.get_acute_status(agent, time)

        if HIVstatus_Agent == 1 and HIVstatus_Partner == 1:
            print "BOTH AGENTS %d and %d WERE HIV+!?!?!?" %(agent._ID, partner._ID)
            return
            #exit(10)
        elif HIVstatus_Agent == 1 or HIVstatus_Partner == 1:
            # Sex between men?
            if Type_agent == 'MSM' and Type_partner == 'MSM':
                SexBetweenMen = 1
            else:
                SexBetweenMen = 0
            # Define probabilities for unsafe sex

            # unprotected sex probabilities for primary partnerships

            # Take average of agent and partner probability
            p_UnsafeSafeSex1 = 0.72 #0.5 * (p_UnsafeSex_Agent1 + p_UnsafeSex_Partner1) ########## OLD METHOD MUST UPDATE TO MATCH RACIAL DISPARITY

            if SexBetweenMen:
                MSexActs = self.MEAN_SMSM_ACTS
                p_UnsafeSafeSex1 = 0.43
            else:
                MSexActs = self.MEAN_S_ACTS

            #print "%s\t%s\t"%(Race_Agent,Type_agent)
            p_UnsafeSafeSex1 = self.ProbTables[Race_Agent][Type_agent]['UNSAFESEX']
            MSexActs = self.ProbTables[Race_Agent][Type_agent]['NUMSexActs'] * self.SEXSCALINGPARAM
            #print "Unsafe:%.5lf\tMSexActs:%.2lf\tOLDMSexActs:%.2lf"%(p_UnsafeSafeSex1,MSexActs,self.MEAN_S_ACTS)
            #print "MSEX",MSexActs
            T_sex_acts1 = int(poisson.rvs(MSexActs, size=1))

            num_int = rel._total_sex_acts
            #Get condom usage
            if num_int < 10:
                if num_int <= 3:
                    p_UnsafeSafeSex1 = 0.603
                else:
                    p_UnsafeSafeSex1 = 0.660
            else:
                p_UnsafeSafeSex1 = 0.826

            #Reduction of risk acts between partners for condom usage
            U_sex_acts1 = T_sex_acts1
            for n in range(U_sex_acts1):
                if random.random() > p_UnsafeSafeSex1:
                    U_sex_acts1 -= 1

            #Reduction of risk acts between partners for PrEP adherence
            U_sex_acts2 = U_sex_acts1
            for n in range(U_sex_acts1):
                if agent._PrEP_bool or partner._PrEP_bool:
                    if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                        if random.random() < 0.96:
                            U_sex_acts2 -= 1
                    else:
                        if random.random() < 0.76:
                            U_sex_acts2 -= 1


            #U_sex_acts1 = T_sex_acts1#int(T_sex_acts1[0] + .5)
            #print "MeanS_act:%.2lf\tT_sex_acts1:%.2lf\tp_UnsafeSex1:%.2lf\tU_sex_acts1:%.2lf"%(MSexActs, T_sex_acts1,p_UnsafeSafeSex1,U_sex_acts1)
            if U_sex_acts2 >= 1:
                # if agent HIV+
                rel._total_sex_acts += U_sex_acts2
                if HIVstatus_Agent == 1 or HIVstatus_Partner == 1:
                    p = self.get_transmission_probability(agent, 'SEX')
                    #p_transmission = binom.pmf(1, U_sex_acts1, p)

                    p_total_transmission = 0
                    if U_sex_acts2 == 1:
                        p_total_transmission = p
                    else:
                        for k in range(1, U_sex_acts2+1):
                            temp = binom.pmf(k, U_sex_acts1, p)
                            #print temp
                            p_total_transmission += temp

                    #print "\t\t\tHIV+ SEX Act A:%d on P:%d\t Must be less than %.10lf\t(k:1 n:%.2lf, p:%.5lf) **OLD pTrans:%.5lf\tAcute:%s" % (agent, partner, p_total_transmission, U_sex_acts1, p, p_transmission, isAcute)
                    #print "\t\t\tHIV+ SEX ACT\tMust be less than %.10lf" % p_total_transmission
                    if random.random() < p_total_transmission:


                        # if agent HIV+ partner becomes HIV+
                        self.Transmit_from_agents += [agent]
                        self.Transmit_to_agents += [partner]
                        if Type_agent == 'MSM': self.Transmission_tracker['SEX_MSM'][time] += 1
                        if Type_agent != 'MSM': self.Transmission_tracker['SEX_NMSM'][time] += 1

                        #print "\t\t\t\tST_HIV from agent %d to partner %d @ %.5lf" %(agent, partner, p)
                        self._become_HIV(partner, time)
                        #print 'INFECTION', Type_agent, DrugType_Agent
                        if(isAcute):
                            self.Acute_agents += [agent]
                            print "\t\t(ACUTE)\tST_HIV from agent %d to partner %d \t@ p=%.5lf transmissionp=%.5lf n:%d" %(agent._ID, partner._ID, p, p_total_transmission, U_sex_acts1)
                        else:
                            print "\t\t\t\tST_HIV from agent %d to partner %d \t@ p=%.5lf transmissionp=%.5lf n:%d" %(agent._ID, partner._ID, p, p_total_transmission, U_sex_acts1)


                        """if agent in ((self.HIV_key_transitiontime[time - 1] if time > 1 else [])
                                         or (self.HIV_key_transitiontime[time - 2] if time > 2 else [])
                                         or (self.HIV_key_transitiontime[time - 3] if time > 3 else [])):
                            self.Acute_agents += [agent]
                            print '\t\t\t\tACUTE', agent"""

                # if partner HIV+
                """if HIVstatus_Agent == 0 and HIVstatus_Partner == 1:
                    p = self.get_transmission_probability(partner, 'SEX')
                    p_transmission = binom.pmf(1, U_sex_acts1, p)
                    #print "\t\t\tHIV+ SEX ACT\tMust be less than %.10lf" % p_transmission
                    if random.random() < p_transmission:
                        # if agent HIV+ partner becomes HIV+
                        self.Transmit_from_agents += [partner]
                        if Type_agent == 'MSM': self.Transmission_tracker['SEX_MSM'][time] += 1
                        if Type_agent != 'MSM': self.Transmission_tracker['SEX_NMSM'][time] += 1

                        print "\t\t\t\tST_HIV from partner %d to agent %d @ %.5lf" %(partner, agent, p)
                        self._become_HIV(agent, time)
                        if partner in ((self.HIV_key_transitiontime[time - 1] if time > 1 else [])
                                           or (self.HIV_key_transitiontime[time - 2] if time > 2 else [])
                                           or (self.HIV_key_transitiontime[time - 3] if time > 3 else [])):
                            self.Acute_agents += [partner]
                            print '\t\t\t\tACUTE', partner"""
            else:
                return
                #print "U_sex = %d < 1" % U_sex_acts1


    #@profile
    def _become_HIV(self, agent, time):
        """
        :Purpose:
            agent becomes HIV agent. Update all appropriate list and
            dictionaries.

        :Input:
            agent : int

        """
        #print "\t!\tAgent %d just became HIV" % agent
        #agent_cl = self.totalAgentClass.get_agent(agent)

        agent._HIV_bool = True
        agent._HIV_time = 1
        self.NewInfections.add_agent(agent)
        print "\t\t\t\tAgent %d added to new infection list"%agent.get_ID()

        self.HIV_agents_class.add_agent(agent)
        #self.totalAgentClass._subset["HIV"].add_agent(agent)
        """
        self.tmp_Agents[agent].update({'HIV': 1})
        if agent not in self.tmp_HIV_agents:
            self.tmp_HIV_agents.append(agent)

        # Keep track of HIV trasmission time
        if self.Agents[agent]['HIV'] == 0:
            self.HIV_key_transitiontime[time].append(agent)

        # write to dynnetworkReport
        #if self.Agents[agent]['HIV'] == 0 or time == 0:
        #    repTime = repr(time) if time > 0 else '1'
        #    reportLine = '\t'.join([repTime, 'BECOMEHIV', repr(agent)])
        #    open('Results/dynnetworkReport.txt', 'a').write('\n' + reportLine)

        #Assign new HIV agent partners
        agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
        agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
        target_num = self._get_number_of_partners(agent, agent_drug_type, agent_sex_type)
        if time > 0:
            AvailableAgents = list(set(self.Agents).difference(set(self.Incarcerated)))
            for n in range(target_num):

                partner = self._get_partner(agent, AvailableAgents)

                if partner != None:
                    partner_sex_type = self.get_agent_characteristic(partner, 'Sex Type')
                    partner_drug_type = self.get_agent_characteristic(partner, 'Drug Type')
                    #print "MATCHING NEW HIV AGENT A:%d (%s %s) with P:%d (%s %s)"%(agent, agent_sex_type, agent_drug_type, partner, partner_sex_type, partner_drug_type)
                    self.AdjMat[agent, partner] = 1
                    self.AdjMat[partner, agent] = 1
                else:
                    print "** COULDNT MATCH AGENT A:%d (%s %s) with anyone"%(agent, agent_sex_type, agent_drug_type)
                    pass

        # keep track of viral load
        adherenceStat = self.AdherenceAgents[agent]
        agent_aids = self.get_agent_characteristic(agent, 'AIDS')
        agent_HIV = 1  # agent_HIV=self.get_agent_characteristic(agent,'HIV')
        agent_acuteStatus = self.get_acute_status(agent, time)
        if agent_aids:
            if adherenceStat == 1:
                Viral_load = 12
            elif adherenceStat == 2:
                Viral_load = 12
            elif adherenceStat == 3:
                Viral_load = 12
            elif adherenceStat == 4:
                Viral_load = 1.4
            elif adherenceStat == 5:
                Viral_load = 12
            elif adherenceStat == 0:
                Viral_load = 7.0
            else:
                print "BAD ADHERENCE STAT!! Exiting program"
                exit(1)

            #print "AIDS PATIENT"

        elif agent_acuteStatus:
            if adherenceStat == 1:
                Viral_load = 6.9
            elif adherenceStat == 2:
                Viral_load = 5.5
            elif adherenceStat == 3:
                Viral_load = 2.6
            elif adherenceStat == 4:
                Viral_load = 1.4
            elif adherenceStat == 5:
                Viral_load = 0.5
            elif adherenceStat == 0:
                Viral_load = 6.9
            else:
                print "BAD ADHERENCE STAT!! Exiting program"
                exit(1)
        else:
            if adherenceStat == 1:
                Viral_load = 4.5
            elif adherenceStat == 2:
                Viral_load = 4.0
            elif adherenceStat == 3:
                Viral_load = 3.3
            elif adherenceStat == 4:
                Viral_load = 2.7
            elif adherenceStat == 5:
                Viral_load = 0.5
            elif adherenceStat == 0:
                Viral_load = 4.5
            else:
                print "BAD ADHERENCE STAT!! Exiting program"
                exit(1)

        self.Viral_load.update({agent: Viral_load})
        """
        #print "\t\t\t\t\tNew HIV agent %d     VL:%.1lf    time:%d"%(agent, Viral_load, time)

    #@profile
    def _update_partner_assignments(self, partnerTurnover):
        # Generate target partner numbers for each agent and get current partner nums
        target_partner_nums = {}
        current_partner_nums = {}

        print "ZZZZZZZZZZZZ"
        print "ZZZZZZZZZZZZ"
        print "ZZZZZZZZZZZZ"
        print "ZZZZZZZZZZZZ"
        print "ZZZZZZZZZZZZ"
        print "ZZZZZZZZZZZZ"

        # TODO: FIX THIS BACK TO HIV AGENTS ONLY
        EligibleAgents = self.Agents#list(set(self.HIV_agents).difference(set(self.Incarcerated)))
        for agent in EligibleAgents:#self.HIV_agents:#Agents:

            agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
            agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')

            current_num = len(list(self.AdjMat.rows[agent]))

            if np.random.uniform(0, 1) > partnerTurnover:
                target_num = current_num
            else:
                target_num = get_number_of_partners(self, agent, agent_drug_type, agent_sex_type)

            target_partner_nums.update({agent: target_num})
            current_partner_nums.update({agent: current_num})

        # Now loop through agents, if currently too many partners, remove some
        for agent in EligibleAgents:#self.HIV_agents:#Agents:
            if target_partner_nums[agent] < current_partner_nums[agent]:
                ExistingLinks = list(self.AdjMat.rows[agent])
                n = current_partner_nums[agent] - target_partner_nums[agent]
                for i in range(n):
                    agent2remove = random.choice(ExistingLinks)
                    #print "Current agent %d has %d partners and wants %d - " %(agent, current_partner_nums[agent],target_partner_nums[agent]), list(self.AdjMat.rows[agent]), "removing %d"%agent2remove
                    self.AdjMat[agent, agent2remove] = 0  # remove connection adjMat
                    self.AdjMat[agent2remove, agent] = 0
                    current_partner_nums[agent] -= 1
                    if agent2remove in EligibleAgents:#self.HIV_agents:
                        current_partner_nums[agent2remove] -= 1

        # Loop through agents again, if too few: go into need_partners set
        need_new_partners = []
        for agent in EligibleAgents:#self.HIV_agents:#Agents:
            if target_partner_nums[agent] > current_partner_nums[agent]:
                need_new_partners.append(agent)

        need_new_partners = list(np.random.permutation(need_new_partners))

        # Now create partnerships until available partnerships are out
        last_list_size = len(need_new_partners)
        iters_at_one_size = 0
        print "\t\t-FINDING MATCHES FOR",len(need_new_partners),"AGENTS IN NEED \t---"
        while len(need_new_partners) > 0:
            #print len(need_new_partners)
            agent = random.choice(need_new_partners)
            agent_cl = self.totalAgentClass.get_agent(agent)
            AvailableAgents = self.Agents
            if self.Incarcerated != []:
                #AvailableAgents.remove(self.Incarcerated)
                AvailableAgents = list(set(self.Agents).difference(set(self.Incarcerated)))
                #print "REMOVED %d from Avialable Lists"%len(self.Incarcerated)

            #for n in AvailableAgents:
            #    print "Agent %d"%n,AvailableAgents[n]
            partner = self._get_partner(agent, AvailableAgents)#self.Agents)
            #partner = self._get_partner(agent, need_new_partners)

            if partner != None:
                partner_cl = self.totalAgentClass.get_agent(partner)
                agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
                agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
                partner_sex_type = self.get_agent_characteristic(partner, 'Sex Type')
                partner_drug_type = self.get_agent_characteristic(partner, 'Drug Type')
                #print "MATCHING A:%d (%s %s) with P:%d (%s %s)"%(agent, agent_sex_type, agent_drug_type, partner, partner_sex_type, partner_drug_type)
                self.AdjMat[agent, partner] = 1
                self.AdjMat[partner, agent] = 1
                current_partner_nums[agent] += 1

                agent_cl.bond(partner_cl)


                if current_partner_nums[agent] == target_partner_nums[agent]:
                    need_new_partners.remove(agent)

                if partner in self.HIV_agents:
                    current_partner_nums[partner] += 1

                    if current_partner_nums[partner] == target_partner_nums[partner]:
                        need_new_partners.remove(partner)
            #else:
                #print "NO MATCH FOUND FOR A:", agent, "( of total", len(need_new_partners),")"
            if len(need_new_partners) == last_list_size:
                iters_at_one_size += 1
            else:
                iters_at_one_size = 0


            if iters_at_one_size > 100: break
            last_list_size = len(need_new_partners)
            if last_list_size == 0: break
        # The remaining partnerless people can remain partnerless :)
        print "\t\t-COULDNT MATCH",len(need_new_partners),"AGENTS IN NEED \t---"
        # Now assign partnerships in the high risk group
        """for agent in self.Agents:
            if agent in self.high_risk_agents and np.random.uniform(0, 1) > self.PARTNERTURNOVER:
                for partner in list(self.AdjMat.rows[agent]):
                    if partner in self.high_risk_agents:
                        self.AdjMat[agent, partner] = 0
                        self.AdjMat[partner, agent] = 0
                agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
                agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
                n_new_partners = self._get_number_of_partners(agent, agent_drug_type, agent_sex_type)
                n_new_partners = int((self.HIGH_RISK_PARTNER_ENHANCEMENT - 1) / 2 * n_new_partners)
            for i in range(n_new_partners):
                partner = self._get_partner(agent, list(self.high_risk_agents))
                if partner != None:
                    self.AdjMat[agent, partner] = 1
                    self.AdjMat[partner, agent] = 1"""


    def _get_number_of_partners(self, agent, agent_drug_type, agent_sex_type):
        """
        :Purpose:
            Get number of partners for a agent.
            Drawn from Poisson distribution.

        :Input:
            agent_drug_type : str
            Either 'IDU', 'NIDU', 'ND'

            agent_sex_type : str
            Either 'HM', 'MSM', 'HF', 'WSW'

        :Output:
            NumPartners : int
            Zero partners possible.
        """
        # Check input
        # Drug type
        if agent_drug_type not in ['IDU', 'NIDU', 'ND']:
            raise ValueError("Invalid drug type! %s" % str(agent_drug_type))
        # Sex type
        if agent_sex_type not in ['HM', 'HF', 'MSM', 'WSW']:
            raise ValueError("Invalid sex type! %s" % str(agent_sex_type))

        agent_race_type = self.get_agent_characteristic(agent, 'Race')

        n_trials = self.ProbTables[agent_race_type][agent_sex_type]['NUMPartn']#5
        p_success = .8

        ##Random number of contacts using negative binomial
        if agent_sex_type == 'WSW':
            # n_trials = 1
            # p_success = 0.8
            RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
        elif agent_sex_type == 'MSM' and agent_drug_type != 'NIDU':
            # n_trials = 1
            # p_success = 0.8
            RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
        elif agent_sex_type == 'MSM' and agent_drug_type == 'NIDU':
            # n_trials = 1
            # p_success = 0.8
            RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
            RandNumCont = int(RandNumCont * 2)
        elif agent_drug_type == 'NIDU':
            # n_trials = 1
            # p_success = 0.8
            RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
        elif agent_drug_type == 'IDU':
            n_trials = 7
            p_success = 0.7
            RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
        elif agent_drug_type == 'ND':
            # n_trials = 1
            # p_success= 0.8
            RandNumCont = np.random.negative_binomial(n_trials, p_success, 1)[0]
        if RandNumCont < 0:
            raise ValueError("Invalid number of contacts!%s" % str(RandNumCont))

        if RandNumCont == 0 and np.random.uniform() < .5:
            RandNumCont = 1
        MEAN_PARTNER_YEAR = self.ProbTables[agent_race_type][agent_sex_type]['NUMPartn']
        RandNumCont = poisson.rvs(MEAN_PARTNER_YEAR, size=1)

        if agent in self.IDU_agents:
            RandNumCont = RandNumCont * 1
        #print "Agent %s\t%s\t%s\tPARTNERS:%d"%(agent_race_type, agent_sex_type, agent_drug_type, RandNumCont)
        RandNumCont = 1 ######################## TEMP FIXER
        return RandNumCont

    #@profile
    def _get_partner(self, agent, need_new_partners):
        """
        :Purpose:
            Get partner for agent.

        :Input:
            agent : int

            need_new_partners: list of available partners

        :Output:
            partner: new partner
        """
        def partner_choice(x):
            intersection = list(set(need_new_partners).intersection(set(x))) #?need_new_partners[:100]?
            #intersection.remove(agent) ### TESTING REMOVING SELF

            if intersection == []: return None
            else: return random.choice(intersection)

        agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
        agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
        ExistingLinks = list(self.AdjMat.rows[agent])
        RandomPartner = None

        #print("Finding partner for agent", agent, agent_sex_type, agent_drug_type)
        if agent_drug_type == 'IDU':
            if random.random() < 0.8:
                # choose from IDU agents
                #print "\tChecking for IDU partner for %d" % agent
                AvailableAgents = list(set(self.IDU_agents).difference(set(self.Incarcerated)))
                RandomPartner = partner_choice(AvailableAgents)
                #print "\tReturned: %s" % RandomPartner
            else:
                RandomPartner = self._get_random_sex_partner(agent, need_new_partners)
        elif agent_drug_type in ('ND','NIDU'):
            RandomPartner = self._get_random_sex_partner(agent, need_new_partners)
        else:
            raise ValueError("Check method _get_partners(). Agent not caught!")
        #print RandomPartner
        if RandomPartner == agent: return None
        else: return RandomPartner


        """if agent_sex_type == 'MSM' and agent_drug_type in ('ND','NIDU'):
            # Either new partner is from MSM only (90% of MSM agents)
            # or partner is chosen from population randomly
            if agent in self.MSM_RandomChoice:
                # Only consider partners which are sex compatible
                RandomPartner = self._get_random_sex_partner(agent, need_new_partners)
            else:
                RandomPartner = partner_choice(self.MSM_agents)

        elif agent_sex_type == 'MSM' and agent_drug_type =='IDU':
            # Either new partner is from MSM only (90% of MSM agents)
            # or partner is chosen from population randomly
            if agent in self.MSM_RandomChoice2:
                # Only consider partners which are sex compatible
                RandomPartner = self._get_random_sex_partner(agent, need_new_partners)
            else:
                RandomPartner = partner_choice(self.MSM_agents)

        elif agent_drug_type == 'IDU' and agent_sex_type !='MSM':
            if random.random() < 0.8:
                # choose from IDU agents
                #print "\tChecking for IDU partner for %d" % agent
                RandomPartner = partner_choice(self.IDU_agents)
                #print "\tReturned: %s" % RandomPartner
            else:
                RandomPartner = self._get_random_sex_partner(agent, need_new_partners)

        elif agent_drug_type == 'NIDU': #incorporate assortative mixing for NIDU
            if random.random() < 0.18:
                # choose from IDU agents
                RandomPartner = partner_choice(self.IDU_agents)
            if random.random() < 0.18 + 0.42:# draw from NIDU
                RandomPartner = partner_choice(self.NIDU_agents)
            else:
                RandomPartner = self._get_random_sex_partner(agent, need_new_partners)

        elif agent_sex_type == 'WSW':
            # Either new partner is from WSW only (50% of WSW agents)
            # or partner is chosen from population randomly
            if agent in self.WSW_RandomChoice:
                # Only consider partners which are sex compatible
                RandomPartner = self._get_random_sex_partner(agent, need_new_partners)
            else:
                RandomPartner = partner_choice(self.WSW_agents)

        elif agent_drug_type =='ND':
            RandomPartner = self._get_random_sex_partner(agent, need_new_partners)
        else:
            raise ValueError("Check method _get_partners(). Agent not caught!")
        if RandomPartner == agent: return None
        else: return RandomPartner"""

    #@profile
    def _get_random_sex_partner(self, agent, need_new_partners):
        """
        :Purpose:
            Get a random partner which is sex compatible

        :Input:
            agent: int
            need_new_partners: list of available partners

        :Output:
            partner : int

        """
        #@profile
        def partner_choice(x):
            intersection = list(set(need_new_partners).intersection(set(x)))
            agent_race_type = self.get_agent_characteristic(agent, 'Race')
            #print agent_race_type
            if agent_race_type == 'WHITE':
                Assortive_intersection = list(set(self.White_agents).intersection(intersection))
                if Assortive_intersection == []: print "Couldnt assortive mix (W), picking suitable agent"
                else: return random.choice(Assortive_intersection)
            elif agent_race_type == 'BLACK':
                Assortive_intersection = list(set(self.Black_agents).intersection(intersection))
                if Assortive_intersection == []:
                    print "Couldnt assortive mix (B), picking suitable agent"
                else:
                    #print Assortive_intersection
                    return random.choice(Assortive_intersection)
            if intersection == []: return None
            else: print "NO PATNAS"#return random.choice(intersection)

        #agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
        agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
        #print "\tChecking for sex partner for %d" % agent

        if agent_sex_type not in ['HM','HF','MSM','WSW']:
            raise ValueError("Invalid sex type! %s"%str(agent_sex_type))
        elif agent_sex_type == 'MSM':
            rv = random.random()
            if rv < 0.91:
                RandomPartner = partner_choice(self.MSM_agents)   # MSM agent
            else:
                RandomPartner = partner_choice(self.HF_agents)  # HF agent
        elif agent_sex_type == 'HM':
            rv = random.random()
            if rv < 1:#0.96:
                RandomPartner = partner_choice(self.HF_agents)   # HF agent
            else:
                RandomPartner = partner_choice(self.WSW_agents)  # WSW agent
        elif agent_sex_type == 'HF':
            rv = random.random()
            if rv < 1:#0.95:
                RandomPartner = partner_choice(self.HM_agents)   # HM agent
            else:
                RandomPartner = partner_choice(self.MSM_agents)  # MSM agent
        elif agent_sex_type == 'WSW':
            rv = random.random()
            if rv < 0.91:
                RandomPartner = partner_choice(self.WSW_agents)   # HM agent
            elif rv < 0.91+0.05:
                RandomPartner = partner_choice(self.MSM_agents)  # MSM agent
            else:
                RandomPartner = partner_choice(self.HM_agents)  # WSW agent
        else:
            raise ValueError("Invalid sex type! %s"%str(agent_sex_type))

        #print "\tReturned: %s" % RandomPartner
        return RandomPartner


    def _sex_possible(self, agent_sex_type, partner_sex_type):
        """
        :Purpose:
        Determine if sex is possible.

        :Input:
        agent_sex_type : str

        partner_sex_type : str

        :Output:
        SexPossible : bool
        """

        # Check input
        if agent_sex_type not in ['HM', 'HF', 'MSM', 'WSW']:
            raise ValueError("Invalid agent_sex_type! %s" % str(agent_sex_type))
        if partner_sex_type not in ['HM', 'HF', 'MSM', 'WSW']:
            raise ValueError("Invalid partner_sex_type! %s" % str(
                partner_sex_type))

        # Sex possible
        if agent_sex_type == 'HM' and partner_sex_type in ['HF', 'WSW']:
            SexPossible = True
        #elif partner_sex_type == 'HM' and agent_sex_type in ['HF', 'WSW']:
        #    SexPossible = True
        elif agent_sex_type == 'MSM' and partner_sex_type in ['MSM', 'WSW', 'HF']:
            SexPossible = True
        #elif partner_sex_type == 'MSM' and agent_sex_type in ['MSM', 'WSW', 'HF']:
        #    SexPossible = True
        elif agent_sex_type == 'WSW' and partner_sex_type in ['MSM', 'WSW', 'HM']:
            SexPossible = True
        #elif partner_sex_type == 'WSW' and agent_sex_type in ['MSM', 'WSW', 'HM']:
        #    SexPossible = True
        elif agent_sex_type == 'HF' and partner_sex_type in ['HM', 'MSM']:
            SexPossible = True
        else:
            SexPossible = False

        if agent_sex_type == 'HM' and partner_sex_type == 'HM' and SexPossible:
            raise ValueError("Check _sex_possible method!")

        return SexPossible


    def _drug_cessation(self, agent, agent_drug_type):
        """
        :Purpose:
            Account for drug cessation of IDU to NIDU and NIDU to ND.

        :Input:
            agent : int

        """
        if agent_drug_type == 'IDU':
            if random.random() < 0.017 / 12 / 2:
                self.tmp_Agents[agent].update({'Drug Type': 'NIDU'})
                if agent in self.tmp_IDU_agents:  # agent might have transitioned already
                    self.tmp_IDU_agents.remove(agent)
                if agent not in self.tmp_NIDU_agents:
                    self.tmp_NIDU_agents.append(agent)
        elif agent_drug_type == 'NIDU':
            if random.random() < 0.017 / 12:
                self.tmp_Agents[agent].update({'Drug Type': 'ND'})
                if agent in self.tmp_NIDU_agents:  # agent might have transitioned already
                    self.tmp_NIDU_agents.remove(agent)
                if agent not in self.tmp_ND_agents:
                    self.tmp_ND_agents.append(agent)
                if agent in self.DrugTreatmentAgents_current:
                    self._exit_drug_treatment(agent)
        else:
            raise ValueError('Drug cessation only valid for IDU and NIDU!')

    def _incarcerate(self, agent, time):
        """
        :Purpose:
            Account for drug cessation of IDU to NIDU and NIDU to ND.

        :Input:
            agent : int

        """
        #agent_dict = self.Agents[agent]
        drug_type = agent._DU#self.get_agent_characteristic(agent, 'Drug Type')
        sex_type = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        race_type = agent._race#self.get_agent_characteristic(agent, 'Race')
        hiv_bool = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        tested = agent._tested#self.get_agent_characteristic(agent, 'Tested')
        incar_t = agent._incar_time#agent_dict['incar_t']
        incar_bool = agent._incar_bool
        haart_bool = agent._HAART_bool

        if incar_bool:#agent in self.Incarcerated:
            #incar_t = agent_dict['incar_t']
            agent._incar_time -= 1
            #self.tmp_Agents[agent].update({'incar_t': incar_t})

            #print "Agent %d has %d months left"%(agent, incar_t)

            #get out if t=0
            if incar_t == 1: #FREE AGENT
                self.IncarceratedClass.remove_agent(agent)
                agent._incar_bool = False
                """except:
                    #print self.Incarcerated
                    self.IncarceratedClass.print_agents()
                    print "AGENT %d was attmpeted to be removed at %d"%(agent.get_ID(),incar_t)
                    exit(9)"""
                #partners = list(self.AdjMat.rows[agent])
                #for partner in enumerate(partners):
                #    self.AdjMat[partner, agent] = 1  # force connection

                agent._mean_num_partners = 12
                agent._highrisk_bool = True
                agent._everhighrisk_bool = True
                agent._highrisk_time = 6
                self.HighriskClass.add_agent(agent)

                if hiv_bool:
                    if haart_bool:
                        if random.random() < 0.43:
                            agent._HAART_bool = True

                        ### FORCE ALL HIGH HAARTS ####
                        #self.AdherenceAgents[agent] = 5
                    else:
                        agent._HAART_bool = False

                        """self.tmp_Agents[agent].update({'HAARTa': 1})
                        self.tmp_HAART_agents.append(agent)

                        self.AdherenceAgents[agent] = 1"""
                        ### END FORCE ####

        elif random.random() < self.ProbTables[race_type][sex_type]['INCAR']*2.5*1.8: ############## SCALED BY 2.5

            toss = random.choice( (1, 2) )
            if toss == 1: #JAIL
                timestay = int(random.triangular(8, 21, 15))
                if hiv_bool and not tested:
                    if random.random() < 0.3:
                        agent._tested = True #self.tmp_Agents[agent].update({'Tested': 1})
                        #self.HIVidentified_agents.append(agent)

            else: #PRISON
                timestay = int(random.triangular(45, 96, 60))
                if hiv_bool and not tested:
                    if random.random() < 1:
                        agent._tested = True #self.tmp_Agents[agent].update({'Tested': 1})
                        #self.HIVidentified_agents.append(agent)

            agent._incar_bool = True
            agent._incar_time = timestay
            self.IncarceratedClass.add_agent(agent)

            #PUT PARTNERS IN HIGH RISK
            for tmpA in agent._partners:
                tmpA._mean_num_partners = 3.1
                tmpA._highrisk_bool = True
                tmpA._everhighrisk_bool = True
                tmpA._highrisk_time = 6
                self.HighriskClass.add_agent(tmpA)
            #print "IMPRISON AGENT %d for %d"%(agent,timestay)
            #self.tmp_Agents[agent].update({'Drug Type': 'ND'})
            #if agent not in self.Incarcerated:
            #self.Incarcerated.append(agent)

            #partners = list(self.AdjMat.rows[agent])
            #if partners:
                #print "Agent %d has partners:" % agent, partners
                #intersection = list(set(need_new_partners[:100]).intersection(set(x)))
                #HIV_partners = list(set(self.AdjMat.rows[agent])).intersection(set(self.tmp_HIV_agents))
                #print HIV_partners
            #    for partner in enumerate(partners):
            #        self.AdjMat[partner, agent] = 0  # remove connection adjMat


    def _drugTest(self, agent, time):
        """
        :Purpose:
            Test the agent for HIV. If detected, add to identified list.

        :Input:
            agent : int
            partner : int
            time : int

        :Output:
            none
        """
        # Drug Treatment
        # SEP

        drug_type = agent._DU#self.get_agent_characteristic(agent, 'Drug Type')
        sex_type = agent._SO#self.get_agent_characteristic(agent, 'Sex Type')
        race_type = agent._race#self.get_agent_characteristic(agent, 'Race')
        hiv_Status = agent._HIV_bool#self.get_agent_characteristic(agent, 'HIV')
        tested = agent._tested#self.get_agent_characteristic(agent, 'Tested')
        #print "Testing agent %d\tHIV:%d\tTested:%d" %(agent, hiv_Status,tested)
        if not tested:
            #print "Testing agent %d\t Needs < %.3lf" %(agent, self.ProbTables[race_type][sex_type]['HIVTEST'])
            test_prob = self.ProbTables[race_type][sex_type]['HIVTEST']
            test_prob = test_prob * 1.2 * self.sensTestFreq# / 5
            if random.random() < test_prob: ###WAS / 10
                #print "\t NEW +"
                agent._tested = True
                #self.tmp_Agents[agent].update({'Tested': 1})
                #self.HIVidentified_agents.append(agent)





    def _VCT(self, agent, time):
        """
        :Purpose:
            Account for voluntary Counseling and Testing(VCT)

        :Input:
            agent : int
            partner : int
            time : int

        :Output:
            none
        """
        # Drug Treatment
        # SEP
        drug_type = self.get_agent_characteristic(agent, 'Drug Type')
        SEPstat = False
        if agent in self.SEPAgents:
            if time == self.SEPAgents[agent]:
                SEPstat = True
        if drug_type == 'IDU':
            if SEPstat:
                if random.random() < self.VCT_NSP:  # !!!!!!!!!!!!!!!!!!!!
                    self.VCTAgents.update({agent: time})
            else:
                if random.random() < self.VCT_NoNSP_IDU:  # !!!!!!!!!!!!!!!!!!!
                    self.VCTAgents.update({agent: time})
        if drug_type == 'NIDU':
            if random.random() < self.VCT_NoNSP_NIDU:  # !!!!!!!!!!!!!!!!!!
                self.VCTAgents.update({agent: time})
        elif agent in self.MSM_agents and random.random() < self.VCT_NoNSP_MSM:  # !
            self.VCTAgents.update({agent: time})
        else:
            if random.random() < self.VCT_NoNSP_EE:  # !!!!!!!!!!!!!!!!!!!!!!!!!
                self.VCTAgents.update({agent: time})


    def _SEP(self, agent, time):
        """
        :Purpose:
            Account for SEP (Syringe Exchange Program) for IDU agents. \n

        :Input:
            time : int

        :Output:
            SEPstat : bool
        """
        if self.get_agent_characteristic(agent, 'Drug Type') != 'IDU':
            raise ValueError("_SEP only valid for IDU agents! agent: %d %s" %
                             (agent, str(self.get_agent_characteristic(agent,
                                                                       'Drug Type'))))
        # Drug treatment increases likelihood of SEP use
        if agent in self.DrugTreatmentAgents_current and self.IDU_agents:

            if random.random() < self.NSP_SAT:
                self.SEPAgents.update({agent: time})
                SEPstat = True
            else:
                SEPstat = False
        elif random.random() < self.NSP_NoSAT:
            self.SEPAgents.update({agent: time})
            SEPstat = True
        else:
            SEPstat = False

        return SEPstat


    def _exit_drug_treatment(self, agent):
        """
        Agent exits drug treament.
        """
        self.DrugTreatmentAgents_past.update({agent:
                                                  self.DrugTreatmentAgents_current[agent]})
        del self.DrugTreatmentAgents_current[agent]


    def _enter_and_exit_drug_treatment(self, agent, time):
        """
        :Purpose:
            Account for drug treatment for IDU agents. \n
            Entering drug treatment is similar to SEP, drug treatment has
            a functional relationship given as follows:
            P(treatment (t+1) | IDU or NIDU) = P(treatment (agent,t) | IDU) if x < N
            else: 0 (x >= 0)
            where N is the total number of treatment slots available. \n
            An agent who was already in drug treatment and relapsed, has a
            pobability twice as strong to reenter drug treatment
            at a later point.

        :Input:
            agent : int

            time : int

        :Output:
            bool
        """
        agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')

        N_TrSpots_Max = 100000  # max number of treatment spots
        if (agent in self.DrugTreatmentAgents_current and random.random() < self.SAT_disc):
            self._exit_drug_treatment(agent)
        elif self.N_TreatmentSpots < N_TrSpots_Max:
            if agent in self.SEPAgents:
                prob = self.SAT_NSP
            elif agent in self.DrugTreatmentAgents_past:
                prob = 0.18
            else:
                if agent_drug_type == 'IDU':
                    prob = self.SAT_NoNSP_IDU
                elif agent_drug_type == 'NIDU':
                    prob = self.SAT_NIDU
                else:
                    mssg = 'Drug treatment only valid for NIDU and IDUs! %s'
                    raise ValueError(mssg % agent_drug_type)
            if random.random() < prob and agent not in self.tmp_ND_agents:
                self.DrugTreatmentAgents_current.update({agent: time})
                self.N_TreatmentSpots += 1
        elif self.N_TreatmentSpots == N_TrSpots_Max:
            pass
        else:
            mssg = 'Check self.N_TreatmentSpots! Max value = 500! %d'
            raise ValueError(mssg % self.N_TreatmentSpots)


    def _initiate_HAART(self, agent, time):
        """
        :Purpose:
            Account for HIV treatment through highly active antiretroviral therapy (HAART).
            HAART was implemented in 1996, hence, there is treatment only after 1996.
            HIV treatment assumes that the agent knows his HIV+ status.

        :Input:
            time : int

        :Output:
            none
        """

        HAART_coverage = 0.30

        # Check valid input
        if agent not in self.HIV_agents:
            print "HIV_agents: ", sorted(self.HIV_agents)
            print "tmp_HIV_agents: ", sorted(self.tmp_HIV_agents)
            print "Agent[agent]", self.Agents[agent]
            try:
                print "tmp_Agent[agent]", self.tmp_Agents[agent]
            except KeyError:
                pass
            raise ValueError("HAART only valid for HIV agents!agent:%s" %
                             str(agent))

        if time == 0:
            agent_haarta = self.get_agent_characteristic(agent, 'HAARTa')
            if agent_haarta == 1 and self.AdherenceAgents[agent] == 0:
                tmp_rnd = random.random()
                HAART_ADH = self.ProbTables['WHITE']['PWID']['HAARTadh']
                if tmp_rnd < HAART_ADH:
                    adherence = 5
                else:
                    adherence = random.randint(1,4)

                self.AdherenceAgents.update({agent: adherence})
        # Determine probability of HIV treatment
        if time > 0:  # HAART implemented 1996######################################################################
            agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
            agent_haarta = self.get_agent_characteristic(agent, 'HAARTa')
            agent_HIV = self.get_agent_characteristic(agent, 'HIV')



            #OLD WAY OF PUTTING AGENTS ON HAART
            agent_Test_bool = self.Agents[agent]['Tested']
            if agent_Test_bool == 1:
                if agent_drug_type == 'IDU':
                    if agent in self.DrugTreatmentAgents_current:
                        prob = self.HAART_SAT_IDU
                    else:
                        prob = self.HAART_NoSAT_IDU
                elif agent_drug_type == 'NIDU':
                    if agent in self.DrugTreatmentAgents_current:
                        prob = self.HAART_SAT_NIDU
                    else:
                        prob = self.HAART_NoSAT_NIDU
                elif agent_drug_type == 'ND':
                    prob = self.HAART_NoSAT_EE
            else:
                prob = 0.0

            # Go on HAART
            if agent_haarta != 1:
                #print prob
                if random.random() < prob:
                    self.tmp_Agents[agent].update({'HAARTa': 1})
                    self.tmp_HAART_agents.append(agent)

                    tmp_rnd = random.random()
                    HAART_ADH = self.ProbTables['WHITE']['PWID']['HAARTadh']
                    if tmp_rnd < HAART_ADH:
                        adherence = 5
                    else:
                        adherence = random.randint(1,4)

                    self.AdherenceAgents.update({agent: adherence})


            #ADHERENCE FORCE#######################################################
            #adherence = 4
            #print "Agent given adherence of %d" % adherence


            """
            if agent_drug_type == 'IDU':
                numHIV_agents = len(set(self.HIVidentified_agents).intersection(set(self.tmp_IDU_agents)))
                numHAART_agents = len(set(self.tmp_HAART_agents).intersection(set(self.tmp_IDU_agents)))
                HAART_coverage = 0.75
            else:
                numHAART_agents = len(self.tmp_HAART_agents)
                #numHIV_agents = len(self.tmp_HIV_agents)
                numHIV_agents = len(self.HIVidentified_agents)
                HAART_coverage = 0.30


            numHIV_IDU_agents = len(set(self.tmp_HAART_agents).intersection(set(self.tmp_IDU_agents)))
            numHAART_IDU_agents = len(set(self.tmp_HAART_agents).intersection(set(self.tmp_IDU_agents)))
            #print "HIV: %d \tHAART:%d"%(numHIV_agents, numHAART_agents)
            if agent_haarta != 1 and numHAART_agents < (HAART_coverage * numHIV_agents):
                #if random.random() < prob / 2:  ##########################################################
                self.tmp_Agents[agent].update({'HAARTa': 1})
                self.tmp_HAART_agents.append(agent)

            #""" #NEWER WAY

            if agent in self.DrugTreatmentAgents_current:
                if agent_drug_type == 'IDU':
                    p_discont = self.HAARTdis_SAT_IDU  # Discontinue therapy
                    tmp_adh_SAT = self.Adhere_SAT_IDU  # Adherence
                elif agent_drug_type == 'NIDU':
                    p_discont = self.HAARTdis_SAT_NIDU
                    tmp_adh_SAT = self.Adhere_SAT_NIDU
                else:
                    mssg = 'Agent: %d\nDrug Type: %s\nDrug treatment entry: %d\n'
                    print(mssg % (agent, agent_drug_type,
                                  self.DrugTreatmentAgents_current[agent]))
                    raise ValueError('Only IDU and NIDU valid!')
            else:
                if agent_drug_type == 'IDU':
                    p_discont = self.HAARTdis_NoSAT_IDU
                    tmp_adh_SAT = self.Adhere_NoSAT_IDU
                elif agent_drug_type == 'NIDU':
                    p_discont = self.HAARTdis_NoSAT_NIDU
                    tmp_adh_SAT = self.Adhere_NoSAT_NIDU
                elif agent_drug_type == 'ND':
                    p_discont = self.HAARTdis_EE
                    tmp_adh_SAT = self.Adhere_EE

            """
            if agent_haarta == 1:
                adhere_1234 = (1 - tmp_adh_SAT) / 4.0
                tmp_rnd = random.random()
                if tmp_rnd < adhere_1234:
                    adherence = 1
                elif tmp_rnd < 2 * adhere_1234:
                    adherence = 2
                elif tmp_rnd < 3 * adhere_1234:
                    adherence = 3
                elif tmp_rnd < 4 * adhere_1234:
                    adherence = 4
                elif tmp_rnd < 4 * adhere_1234 + tmp_adh_SAT:
                    adherence = 5

                #adherence = 5###########################################
            else:
                adherence = 0

            #ADHERENCE FORCE#######################################################
            #adherence = 4
            #print "Agent given adherence of %d" % adherence
            self.AdherenceAgents.update({agent: adherence})

            """

            if agent_haarta == 1 and random.random() < 0: #p_discont: ###########################################
                self.HAART_agents.remove(agent)
                self.AdherenceAgents.update({agent: 0})

    def _initiate_PrEP(self, agent, time):
        """
        :Purpose:
            Place agents onto PrEP treatment.
            PrEP treatment assumes that the agent knows his HIV+ status is negative.

        :Input:
            time : int

        :Output:
            none
        """

        # Check valid input
        if agent._HIV_bool:
            raise ValueError("PrEP only valid for HIV- agents!agent:%d" % agent.get_ID())

        # Determine probability of HIV treatment
        agent_drug_type = agent._DU



        #OLD WAY OF PUTTING AGENTS ON PrEP
        if not agent._PrEP_bool:
            if agent_drug_type == 'IDU':
                if agent in self.DrugTreatmentAgents_current:
                    prob = self.HAART_SAT_IDU
                else:
                    prob = self.HAART_NoSAT_IDU
            elif agent_drug_type == 'NIDU':
                if agent in self.DrugTreatmentAgents_current:
                    prob = self.HAART_SAT_NIDU
                else:
                    prob = self.HAART_NoSAT_NIDU
            elif agent_drug_type == 'ND':
                prob = self.HAART_NoSAT_EE
        else:
            prob = 0.0

        # Go on PrEP
        prob = 0.0
        #print prob


        if agent_drug_type == 'IDU':
            numHIV_agents = len(set(self.HIVidentified_agents).intersection(set(self.tmp_IDU_agents)))
            numHAART_agents = len(set(self.tmp_HAART_agents).intersection(set(self.tmp_IDU_agents)))
            HAART_coverage = 0.75
        else:
            numPrEP_agents = self.totalAgentClass._subset["PrEP"].num_members()
            #numHIV_agents = len(self.tmp_HIV_agents)
            numHIV_agents = self.HIVidentified_agents
            HAART_coverage = 0.30
            target_PrEP = 71 *(float(time)/12)

            target_PrEP = 0.0973993*time*time + 0.778751 * time - 7.40869
            target_PrEP = self.totalAgentClass.num_members()*0.01*time/12

        #print "HIV: %d \tHAART:%d"%(numHIV_agents, numHAART_agents)
        if numPrEP_agents < target_PrEP:
            #if random.random() < prob / 2:  ##########################################################
            agent._PrEP_bool = True
            self.totalAgentClass._subset["PrEP"].add_agent(agent)
            tmp_rnd = random.random()
            PrEP_Adh = 0.82
            if tmp_rnd < PrEP_Adh:
                agent._PrEP_adh = 1
            else:
                agent._PrEP_adh = 0

        #""" #NEWER WAY




    def _progress_to_AIDS(self, agent, agent_drug_type):
        """
        :Purpose:
            Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        if not agent._HIV_bool:
            raise ValueError("HAART only valid for HIV agents!agent:%s" % str(agent._ID))

        #if agent not in self.AIDS_agents:
        if not agent._HAART_bool:
            adherenceStat = agent._HAART_adh#self.AdherenceAgents[agent]
            if adherenceStat > 0:
                if adherenceStat == 1:
                    if agent_drug_type == 'IDU':
                        prob = self.PAIDS_IDU_1
                    else:
                        prob = self.PAIDS_EE_1

                    prob = 0.0051
                if adherenceStat == 2:
                    if agent_drug_type == 'IDU':
                        prob = self.PAIDS_IDU_2
                    else:
                        prob = self.PAIDS_EE_2

                    prob = 0.0039

                if adherenceStat == 3:
                    if agent_drug_type == 'IDU':
                        prob = self.PAIDS_IDU_3
                    else:
                        prob = self.PAIDS_EE_3

                    prob = 0.0032

                if adherenceStat == 4:
                    if agent_drug_type == 'IDU':
                        prob = self.PAIDS_IDU_4
                    else:
                        prob = self.PAIDS_EE_4

                    prob = 0.0025

                if adherenceStat == 5:
                    if agent_drug_type == 'IDU':
                        prob = self.PAIDS_IDU_5
                    else:
                        prob = self.PAIDS_EE_5

                    prob = 0.0008

            else:
                if agent_drug_type == 'IDU':
                    prob = self.PAIDS_IDU_0
                else:
                    prob = self.PAIDS_EE_0

                prob = 0.0051

            if random.random() < prob * 0.5:#######################################################
                #print "--------AIDS? -> yes"
                agent._AIDS_bool = True
                #self.tmp_AIDS_agents.append(agent)
                #self.tmp_Agents[agent].update({'AIDS': 1})


    def _reset_death_count(self):
        self.num_Deaths = {}
        for HIV_status in ['Total','HIV-', 'HIV+']:
            self.num_Deaths.update({HIV_status: {}})
            for tmp_type in [HIV_status, 'MSM', 'HM', 'HF', 'WSW']:
                self.num_Deaths[HIV_status].update({tmp_type: 0})


    def _remove_agent(self, agent):
        """
        :Purpose:
            Remove agent from the population.
            Delete agent is a key to an associated dictionary which stores the internal.

        :Input:
            agent : int

        """

        # Drug type lists (agent might have been updated)
        drug_type = self.tmp_Agents[agent]['Drug Type']

        if drug_type == 'IDU':
            try:
                self.tmp_IDU_agents.remove(agent)
            except ValueError:
                pass
        elif drug_type == 'NIDU':
            try:
                self.tmp_NIDU_agents.remove(agent)
            except ValueError:
                pass
        elif drug_type == 'ND':
            try:
                self.tmp_ND_agents.remove(agent)
            except ValueError:
                pass
                #print "Agents[agent]", self.Agents[agent]
                #print "tmp_Agents[agent]", self.tmp_Agents[agent]
                #print "tmp_ND_agents", sorted(self.tmp_ND_agents)
        else:
            raise ValueError("Invalid drug type! %s" % str(drug_type))

        # Sex type lists
        try:
            # Agent might have been updated
            sex_type = self.tmp_Agents[agent]['Sex Type']
        except KeyError:
            sex_type = self.get_agent_characteristic(agent, 'Sex Type')
        if sex_type == 'MSM':
            self.tmp_MSM_agents.remove(agent)
        elif sex_type == 'HF':
            self.tmp_HF_agents.remove(agent)
        elif sex_type == 'WSW':
            self.tmp_WSW_agents.remove(agent)
        elif sex_type == 'HM':
            self.tmp_HM_agents.remove(agent)
        else:
            raise ValueError("Invalid sex type! %s" % str(sex_type))

        # HIV and AIDS lists
        try:
            self.tmp_HIV_agents.remove(agent)
        except ValueError:
            pass

        try:
            self.tmp_AIDS_agents.remove(agent)
        except ValueError:
            pass

        try:
            self.Incarcerated.remove(agent)
        except ValueError:
            pass

        try:
            self.tmp_HAART_agents.remove(agent)
        except ValueError:
            #print "WTFFFFF"
            pass

        try:
            self.HIVidentified_agents.remove(agent)
        except ValueError:
            #print "WTFFFFF"
            pass

        # Other lists / dictionaries
        if agent in self.SEPAgents:
            # dict of users who used SEP (agent:time)
            del self.SEPAgents[agent]
        if agent in self.SEPAgents_past:
            del self.SEPAgents_past[agent]
        if agent in self.DrugTreatmentAgents_current:
            # dictionary of users who are currently undergoing
            del self.DrugTreatmentAgents_current[agent]
        if agent in self.DrugTreatmentAgents_past:
            # dictionary of users who underwent drug treatment
            del self.DrugTreatmentAgents_past[agent]
        if agent in self.VCTAgents:
            # list of agents who get tested for HIV ((agent:time)
            del self.VCTAgents[agent]
        for time in self.HIV_key_transitiontime:
            tmp_agents = self.HIV_key_transitiontime[time]
            if agent in tmp_agents:
                tmp_agents.remove(agent)
            self.HIV_key_transitiontime.update({time: tmp_agents})


    def _die_and_replace(self):
        """
        :Purpose:
            Let agents die and replace the dead agent with a new agent randomly.
        """
        totalDeaths = 0

        #dynnetworkReport = open('Results/dynnetworkReport.txt', 'a')
        for agent in self.totalAgentClass.iter_agents(): #self.Agents:

            if agent._incar_bool:#self.IncarceratedClass.is_member(agent):
                #print "Agent %d is incarcerated. Cannot die" % agent.get_ID()
                pass
            else:
                # Probability for dying
                drug_type = agent._DU #self.get_agent_characteristic(agent, 'Drug Type')
                sex_type = agent._SO #self.get_agent_characteristic(agent, 'Sex Type')
                HIV_status = agent._HIV_bool #self.get_agent_characteristic(agent, 'HIV')
                AIDSStatus = agent._AIDS_bool #self.get_agent_characteristic(agent, 'HIV')
                agent_Race = agent._race #self.get_agent_characteristic(agent, 'Race')
                #adherence =  self.AdherenceAgents[agent]


                if HIV_status:
                    if AIDSStatus: #AIDS DEATH RATE
                        if agent_Race == 'WHITE':
                            p = 21.7
                        elif agent_Race == 'BLACK':
                            p = 33.0
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                        #p = self.ProbDeath[drug_type]['AIDS']

                    elif agent._HAART_adh > 100: #HAART DEATH RATE
                        if agent_Race == 'WHITE':
                            p = 8.89
                        elif agent_Race == 'BLACK':
                            p = 8.8
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                        #p = self.ProbDeath[drug_type]['HIV+/HAART']

                    else: #HIV+ DEATH RATE
                        if agent_Race == 'WHITE':
                            p = 6.28
                        elif agent_Race == 'BLACK':
                            p = 2500#32.6#16.5
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                        #p = self.ProbDeath[drug_type]['HIV+']

                elif not HIV_status: # NON HIV DEATH RATE
                    if agent_Race == 'WHITE':
                        p = 3.6#7.9
                    elif agent_Race == 'BLACK':
                        p = 0#8.8
                    else:
                        raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                    #p = self.ProbDeath[drug_type]['HIV-']

                else:
                    raise ValueError("Invalid HIV type! %s" % str(HIV_status))

                #print p
                p = p / 12000.0#12000.0 #putting it into per 1 person-month

                #p=0#****** MAKE SURE TO FIX THIS AFTER!!!! JUST FOR DEBUG PURPOSES!!!!!********
                if random.random() < p:
                    print "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tAgent %d died rolling under %.10lf" % (agent.get_ID(), p)
                    agent.print_agent()
                    totalDeaths += 1
                    if HIV_status: ident = "HIV+"
                    else: ident = "HIV-"
                    self.num_Deaths["Total"][sex_type] += 1
                    self.num_Deaths[ident][sex_type] += 1
                    ID_number = agent.get_ID()
                    race = agent._race

                    #self._remove_agent(agent)
                    #self.totalAgentClass.print_agents()
                    self.totalAgentClass.remove_agent(agent)
                    del agent
                    #agent.print_agent()
                    agent_cl = self._return_new_Agent_class(ID_number,race)
                    agent_cl._HIV_bool = False
                    agent_cl._HIV_time = 0
                    agent_cl._AIDS_bool = False
                    agent_cl._AIDS_time = 0
                    agent_cl._tested = False
                    agent_cl._HAART_bool = False
                    agent_cl._HAART_time = 0
                    agent_cl._HAART_adh = 0
                    agent_cl._PrEP_bool = False
                    agent_cl._PrEP_time = 0
                    agent_cl._PrEP_adh = 0
                    agent_cl._incar_bool = False
                    agent_cl._incar_time = 0
                    #print "adding new agent to totalagents"
                    self.totalAgentClass.add_agent(agent_cl)

                elif 1==0:
                    # Replace with new agent (random characteristics)
                    rv = random.random()
                    if rv < 0.571:
                        deliminator = 'WHITE'
                        #print "\t\tReplaced with ND"
                    else:
                        deliminator = 'BLACK'
                        #print "\t\tReplaced with NIDU"

                        """#################### OLD WAY
                    # Replace with new agent (random characteristics)
                    rv = random.random()
                    if rv < 0.9229:
                        drug_type = 'ND'
                        print "\t\tReplaced with ND"
                    elif rv < 0.9229 + 0.0647:
                        drug_type = 'NIDU'
                        print "\t\tReplaced with NIDU"
                    else:
                        drug_type = 'IDU'
                        print "\t\tReplaced with "
                        """


                     #################### NOW SET TO REPLAC WITH WHAT DIED"

                    #drug_type = 'IDU'

                    # New agent dict
                    agent_dict = self._return_new_agent_dict(deliminator)
                    if deliminator != agent_dict['Race']:
                        raise ValueError("Inconsistent drug type!%s" % str(agent_dict['Drug Type']))

                    # Update tmp_Agents dictionary with new agent
                    self.tmp_Agents.update({agent: agent_dict})
                    #print "New agent updated"

                    # Drug Type
                    drug_type = agent_dict['Drug Type']
                    if drug_type == 'IDU':
                        self.tmp_IDU_agents.append(agent)
                    elif drug_type == 'NIDU':
                        self.tmp_NIDU_agents.append(agent)
                    elif drug_type == 'ND':
                        self.tmp_ND_agents.append(agent)
                    else:
                        raise ValueError("Invalid drug type! %s" % str(drug_type))

                    # Sex Type
                    SexType = agent_dict['Sex Type']
                    if SexType == 'HM':
                        self.tmp_HM_agents.append(agent)
                    elif SexType == 'HF':
                        self.tmp_HF_agents.append(agent)
                    elif SexType == 'MSM':
                        self.tmp_MSM_agents.append(agent)
                    elif SexType == 'WSW':
                        self.tmp_WSW_agents.append(agent)
                    else:
                        raise ValueError("Invalid SexType! %s" % str(SexType))

                    # HIV
                    HIVStatus = agent_dict['HIV']
                    if HIVStatus == 1:
                        #print "NEW AGENT %d %s WAS HIV"%(agent, drug_type)
                        self.tmp_HIV_agents.append(agent)
                    elif HIVStatus != 0:
                        raise ValueError("Invalid HIVType! %s" % str(HIVStatus))
                    #else:
                        #print "NEW AGENT %d %s WAS NOT HIV"%(agent, drug_type)

                    # AIDS
                    AIDSStatus = agent_dict['AIDS']
                    if AIDSStatus == 1:
                        #print "NEW AGENT WAS AIDS"
                        self.tmp_AIDS_agents.append(agent)
                    elif AIDSStatus != 0:
                        raise ValueError("Invalid AIDS Status! %s" % str(AIDSStatus))

                    # HAART
                    HAARTStatus = agent_dict['HAARTa']
                    if HAARTStatus == 1:
                        #print "NEW AGENT WAS HAART"
                        self.tmp_HAART_agents.append(agent)
                    elif HAARTStatus != 0:
                        raise ValueError("Invalid HAART Status! %s" % str(HAARTStatus))

                    #Incarcerated
                    IncarceratedTime = agent_dict['incar_t']
                    if IncarceratedTime >= 1:
                        self.Incarcerated.append(agent)
                    elif IncarceratedTime < 0:
                        raise ValueError("Invalid AIDS Status! %s"%str(IncarceratedTime))

                    # Check
                    if HIVStatus == 1:
                        if agent not in self.tmp_HIV_agents:
                            raise ValueError("Agent must be in HIV_agents")
                    if AIDSStatus == 1:
                        if agent not in self.tmp_AIDS_agents:
                            raise ValueError("Agent must be in AIDS_agents")
                    if HAARTStatus == 1:
                        if agent not in self.tmp_HAART_agents:
                            raise ValueError("Agent must be in HAART_agents")

                    """# write new agent to dynnetworkReport
                    #print "Writing to dynNetReport"
                    reportLine = '\t'.join([repr(self.TimeStep), 'DEATH', repr(agent)])
                    #open('Results/dynnetworkReport.txt', 'a').write('\n' + reportLine)
                    dynnetworkReport.write('\n' + reportLine)

                    reportLine = '\t'.join([repr(self.TimeStep), 'NEWAGENT', repr(agent), SexType, drug_type, repr(HIVStatus)])
                    #open('Results/dynnetworkReport.txt', 'a').write('\n' + reportLine)
                    dynnetworkReport.write('\n' + reportLine)"""

            #if len(self.Agents) != self.PopulationSize:
            #    raise ValueError("Wrong Population size!%s != %s" % (str(len(self.Agents)), str(self.PopulationSize)))

            #iduPrec = self.num_Deaths['IDU']['IDU']/float(totalDeaths) if totalDeaths else 0
            #niduPrec = self.num_Deaths['NIDU']['NIDU']/float(totalDeaths) if totalDeaths else 0
            #ndPrec = self.num_Deaths['ND']['ND']/float(totalDeaths) if totalDeaths else 0
            """
            print "\n\t=== DEATH REPORT ==="
            print "\tType\tCount\t%"
            print "\tIDU \t%d\t\t%.5lf" % (self.num_Deaths['IDU']['IDU'], iduPrec)
            print "\tNIDU\t%d\t\t%.5lf" % (self.num_Deaths['NIDU']['NIDU'], niduPrec)
            print "\tND  \t%d\t\t%.5lf" % (self.num_Deaths['ND']['ND'], ndPrec)
            print "\tTotal\t%d\n" % totalDeaths
            """
            #dynnetworkReport.close()

    #@profile
    def _update_population(self):
        """
        :Purpose:
            Update the population. Changes resulting from parsing
            through the agents and applying the update rules are stored
            in :py:attr:`tmp_agent_dict`. This method updates the whole
            population, i.e., it copies changes from the :py:attr:`tmp_agent_dict`
            dictionary and copies it into the :py:attr:`Agents` dictionary.
        :Input:
            none

        :Output:
            none
        """

        # self.Agents = deepcopy(self.tmp_Agents)
        #
        # self.IDU_agents = copy(list(set(self.tmp_IDU_agents)))
        # self.NIDU_agents = copy(list(set(self.tmp_NIDU_agents)))
        # self.ND_agents = copy(list(set(self.tmp_ND_agents)))
        #
        # self.HM_agents = copy(list(set(self.tmp_HM_agents)))
        # self.HF_agents = copy(list(set(self.tmp_HF_agents)))
        # self.MSM_agents = copy(list(set(self.tmp_MSM_agents)))
        # self.WSW_agents = copy(list(set(self.tmp_WSW_agents)))
        #
        # self.AIDS_agents = copy(list(set(self.tmp_AIDS_agents)))
        # self.HIV_agents = copy(list(set(self.tmp_HIV_agents)))
        # self.HAART_agents = copy(list(set(self.tmp_HAART_agents)))

        #print "HIV population:", self.HIV_agents
        self.SEPAgents = {}  # SEP has no memory


    def _check_population(self):
        """
        :Purpose:
            Check consistency of population.
            Only called in unittest.

        """

        # Check consistency of last partners
        if (not (np.all(self.AdjMat.sum(0) == self.AdjMat.conj().sum(0)) and
                     np.all(self.AdjMat.sum(1) == self.AdjMat.conj().sum(1)))):
            raise ValueError("Adjacency matrix not symmetric!")

        # Check consistency of real population
        count_HF = 0
        count_HM = 0
        count_MSM = 0
        count_WSW = 0
        count_ND = 0
        count_NIDU = 0
        count_IDU = 0
        count_HIV = 0
        count_AIDS = 0
        for (agent, d) in self.Agents.iteritems():
            agent_dict = d
            # Sex type
            sex_type = agent_dict['Sex Type']
            if sex_type == 'HF':
                if agent not in self.HF_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents HF Sex type %d" % agent)
                else:
                    count_HF += 1
            elif sex_type == 'HM':
                if agent not in self.HM_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents HM Sex type %d" % agent)
                else:
                    count_HM += 1
            elif sex_type == 'MSM':
                if agent not in self.MSM_agents:
                    raise ValueError("Check agents MSM Sex type %d" % agent)
                else:
                    count_MSM += 1
            elif sex_type == 'WSW':
                if agent not in self.WSW_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents WSW Sex type %d" % agent)
                else:
                    count_WSW += 1
            else:
                raise ValueError("Invalid sex type %s" % str(sex_type))

            # Drug type
            drug_type = agent_dict['Drug Type']
            if drug_type == 'ND':
                if agent not in self.ND_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents ND Drug type %d" % agent)
                else:
                    count_ND += 1
            elif drug_type == 'NIDU':
                if agent not in self.NIDU_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents NIDU Drug type %d" % agent)
                else:
                    count_NIDU += 1
            elif drug_type == 'IDU':
                if agent not in self.IDU_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agents IDU Drug type %d" % agent)
                else:
                    count_IDU += 1
            else:
                raise ValueError("Invalid drug type %s" % str(drug_type))

            # HIV
            HIVstatus = agent_dict['HIV']
            if HIVstatus != 0:
                if agent not in self.HIV_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agent HIV %d" % agent)
                else:
                    count_HIV += 1
            # AIDS
            AIDSstatus = agent_dict['AIDS']
            if AIDSstatus != 0:
                if agent not in self.AIDS_agents:
                    print self.Agents[agent]
                    raise ValueError("Check agent AIDS %d" % agent)
                else:
                    count_AIDS += 1

        if len(self.HF_agents) != count_HF:
            raise ValueError("self.HF agents contains too many agents!")
        if len(self.HM_agents) != count_HM:
            print "len(self.HM_agents)=%d" % len(self.HM_agents)
            print "count_HM=%d" % count_HM
            raise ValueError("self.HM agents contains too many agents!")
        if len(self.MSM_agents) != count_MSM:
            raise ValueError("self.MSM agents contains too many agents!")
        if len(self.WSW_agents) != count_WSW:
            raise ValueError("self.WSW agents contains too many agents!")

        if len(self.NIDU_agents) != count_NIDU:
            raise ValueError("self.NIDU_agents contains too many agents!")
        if len(self.ND_agents) != count_ND:
            raise ValueError("self.ND agents contains too many agents!")
        if len(self.IDU_agents) != count_IDU:
            mssg = "self.IDU agents contains too many agents!\
                    \nlen(self.IDU_agents)=%d\ncount_IDU=%d\n"
            raise ValueError(mssg % (len(self.IDU_agents), count_IDU))

        if len(self.HIV_agents) != count_HIV:
            raise ValueError("self.HIV_agents contains too many agents!\
                \nlen(self.HIV_agents) = %d\ncount_HIV = %d\n" % (
                len(self.HIV_agents), count_HIV))
        if len(self.AIDS_agents) != count_AIDS:
            raise ValueError("self.AIDS agents contains too many agents!")

        # Check consistency of tmp population
        count_HF = 0
        count_HM = 0
        count_MSM = 0
        count_WSW = 0
        count_ND = 0
        count_NIDU = 0
        count_IDU = 0
        count_HIV = 0
        count_AIDS = 0
        for (agent, d) in self.tmp_Agents.iteritems():
            agent_dict = d
            # Sex type
            sex_type = agent_dict['Sex Type']
            if sex_type == 'HF':
                if agent not in self.tmp_HF_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_HF += 1
            elif sex_type == 'HM':
                if agent not in self.tmp_HM_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_HM += 1
            elif sex_type == 'MSM':
                if agent not in self.tmp_MSM_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_MSM += 1
            elif sex_type == 'WSW':
                if agent not in self.tmp_WSW_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Sex type %d" % agent)
                else:
                    count_WSW += 1
            else:
                raise ValueError("Invalid sex type %s" % str(sex_type))

            # Drug type
            drug_type = agent_dict['Drug Type']
            if drug_type == 'ND':
                if agent not in self.tmp_ND_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Drug type %d" % agent)
                else:
                    count_ND += 1
            elif drug_type == 'NIDU':
                if agent not in self.tmp_NIDU_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Drug type %d" % agent)
                else:
                    count_NIDU += 1
            elif drug_type == 'IDU':
                if agent not in self.tmp_IDU_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agents Drug type %d" % agent)
                else:
                    count_IDU += 1
            else:
                raise ValueError("Invalid drug type %s" % str(drug_type))

            # HIV
            HIVstatus = agent_dict['HIV']
            if HIVstatus != 0:
                if agent not in self.tmp_HIV_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check tmp_agent HIV %d" % agent)
                else:
                    count_HIV += 1
            # AIDS
            AIDSstatus = agent_dict['AIDS']
            if AIDSstatus != 0:
                if agent not in self.tmp_AIDS_agents:
                    print self.tmp_Agents[agent]
                    raise ValueError("Check agent AIDS %d" % agent)
                else:
                    count_AIDS += 1

        if len(self.tmp_HF_agents) != count_HF:
            raise ValueError("self.tmp_HF agents contains too many agents!")
        if len(self.tmp_HM_agents) != count_HM:
            raise ValueError("self.tmp_HM agents contains too many agents!")
        if len(self.tmp_MSM_agents) != count_MSM:
            raise ValueError("self.tmp_MSM agents contains too many agents!")
        if len(self.tmp_WSW_agents) != count_WSW:
            raise ValueError("self.tmp_WSW agents contains too many agents!")

        if len(self.tmp_NIDU_agents) != count_NIDU:
            raise ValueError("self.tmp_NIDU_agents contains too many agents!")
        if len(self.tmp_ND_agents) != count_ND:
            raise ValueError("self.tmp_ND agents contains too many agents!")
        if len(self.tmp_IDU_agents) != count_IDU:
            mssg = "self.tmp_IDU agents contains too many agents!\
                    \nlen(self.tmp_IDU_agents)=%d\ncount_IDU=%d\n"
            raise ValueError(mssg % (len(self.IDU_agents), count_IDU))

        if len(self.tmp_HIV_agents) != count_HIV:
            raise ValueError("self.tmp_HIV_agents contains too many agents!")
        if len(self.tmp_AIDS_agents) != count_AIDS:
            print "len(self.tmp_AIDS_agents)=%d" % len(self.tmp_AIDS_agents)
            print "count_AIDS=%d" % count_AIDS
            raise ValueError("self.tmp_AIDS agents contains too many agents!")


    def save_AgentPartner_list(self, t):
        """
        :Purpsose:
        Save all agent-partners connections.
        :Input:
        t : int
        Time
        """
        OutFileDir = os.path.expanduser(os.path.join(self.current_dir, 'Results'))
        if not os.path.isdir(OutFileDir):  # create directory if not existing
            os.mkdir(OutFileDir)
        OutFileName = os.path.join(OutFileDir,
                                   'AgentPartnersList_atTime_%s.txt' % str(t))
        if os.path.isfile(OutFileName): os.remove(OutFileName)
        outfile = open(OutFileName, 'w')
        outfile.write('agent\tdrug type\tsex type\tHIV\tAIDS\tHAART\t')
        maxpartners = 0
        for agent in self.Agents:
            numpartners = len(list(self.AdjMat.rows[agent]))
            if numpartners > maxpartners:
                maxpartners = numpartners
        outfile.write('\t'.join(['partner\tp drug type\tp sex type'] *
                                maxpartners))
        outfile.write('\n')
        for agent in sorted(self.Agents.keys()):
            agent_dict = self.Agents[agent]
            outfile.write('%d\t' % agent)
            outfile.write('%s\t' % agent_dict['Drug Type'])
            outfile.write('%s\t' % agent_dict['Sex Type'])
            outfile.write('%d\t' % agent_dict['HIV'])
            outfile.write('%d\t' % agent_dict['AIDS'])
            outfile.write('%d\t' % self.AdherenceAgents[agent])
            for p in sorted(list(self.AdjMat.rows[agent])):
                partner_dict = self.Agents[p]
            outfile.write('%d\t' % int(p))
            outfile.write('%s\t' % partner_dict['Drug Type'])
            outfile.write('%s\t' % partner_dict['Sex Type'])
            outfile.write('\n')


    def _reset_partner_count(self):
        """
        Reset partner count for method assess_interaction_distribution
        """

        # set ND partner count to zero for the next time step
        self.tmp_ND_NumPartners_Count = {}
        self.tmp_NIDU_NumPartners_Count = {}
        self.tmp_IDU_NumPartners_Count = {}
        self.tmp_MSM_NumPartners_Count = {}


    def get_HIV_prevalence_drugs(self):
        """
        get HIV prevalence within all three drug user groups
        """
        count_HIV_IDU = 0
        count_HIV_NIDU = 0
        count_HIV_ND = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_characteristic(agent, 'HIV')
            if HIVstatus == 1:
                agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
            if agent_drug_type == 'IDU':
                count_HIV_IDU += 1
            elif agent_drug_type == 'NIDU':
                count_HIV_NIDU += 1
            elif agent_drug_type == 'ND':
                count_HIV_ND += 1
            elif HIVstatus != 0:
                print HIVstatus
                raise ValueError("HIV status must be either 0 or 1 !")
                # print [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]
            else:
                raise ValueError("Agent must be either IDU, NIDU or ND !")
        return [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]


    def get_HIV_prevalence_sex(self):
        """ get HIV prevalence within all four sex groups """
        count_HIV_MSM = 0
        count_HIV_HM = 0
        count_HIV_HF = 0
        count_HIV_WSW = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_characteristic(agent, 'HIV')
            if HIVstatus == 1:
                agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
            if agent_sex_type == 'MSM':
                count_HIV_MSM += 1
            elif agent_sex_type == 'HM':
                count_HIV_HM += 1
            elif agent_sex_type == 'HF':
                count_HIV_HF += 1
            elif agent_sex_type == 'WSW':
                count_HIV_WSW += 1
            elif HIVstatus != 0:
                print HIVstatus
                raise ValueError("HIV status must be either 0 or 1 !")
                # print [count_HIV_IDU, count_HIV_NIDU, count_HIV_ND]
            else:
                raise ValueError("Agent must be either MSM, HM, MF, or WSW !")

        return [count_HIV_MSM, count_HIV_HM, count_HIV_HF, count_HIV_WSW]


    def get_HIV_prevalence_drugs_sex(self):
        """prevalences without and msm only"""
        count_HIV_MIDU = 0
        count_HIV_MNIDU = 0
        count_HIV_MND = 0
        count_HIV_IDUnmsm = 0
        count_HIV_NIDUnmsm = 0
        count_HIV_NDnmsm = 0

        for agent in self.Agents:
            HIVstatus = self.get_agent_characteristic(agent, 'HIV')
            if HIVstatus == 1:
                agent_sex_type = self.get_agent_characteristic(agent, 'Sex Type')
            agent_drug_type = self.get_agent_characteristic(agent, 'Drug Type')
            if agent_drug_type == 'IDU' and agent_sex_type in ['HM', 'HF', 'WSW']:
                count_HIV_IDUnmsm += 1
            elif agent_drug_type == 'IDU' and agent_sex_type == 'MSM':
                count_HIV_MIDU += 1
            elif agent_drug_type == 'NIDU' and agent_sex_type in ['HM', 'HF', 'WSW']:
                count_HIV_NIDUnmsm += 1
            elif agent_drug_type == 'NIDU' and agent_sex_type == 'MSM':
                count_HIV_MNIDU += 1
            elif agent_drug_type == 'ND' and agent_sex_type in ['HM', 'HF', 'WSW']:
                count_HIV_NDnmsm += 1
            elif agent_drug_type == 'ND' and agent_sex_type == 'MSM':
                count_HIV_MND += 1
            elif HIVstatus != 0:
                print HIVstatus
            raise ValueError("HIV status must be either 0 or 1 !")
        return [count_HIV_MIDU, count_HIV_MNIDU, count_HIV_MND, count_HIV_IDUnmsm, count_HIV_NIDUnmsm, count_HIV_NDnmsm]


    def get_HIV_prevalence(self):
        """ get HIV prevalence"""
        HIVcount = 0.0
        for agent in self.Agents.keys():
            HIVstatus = self.get_agent_characteristic(agent, 'HIV')
            if HIVstatus == 1: HIVcount += 1
        return HIVcount


    def return_results(self):
        return self.ResultDict


    def save_result_dict(self):
        OutFileDir = os.path.join(self.current_dir, 'Results')
        if not os.path.isdir(OutFileDir):  # create directory if not existing
            os.mkdir(OutFileDir)
        OutFileName = os.path.join(OutFileDir, 'ResultDictionary.txt')
        if os.path.isfile(OutFileName): os.remove(OutFileName)
        outfile = open(OutFileName, 'w')
        for result_property in sorted(self.ResultDict.keys()):
            outfile.write('%s\t' % result_property)
            for time_t in sorted(self.ResultDict[result_property].keys()):
                outfile.write('%4.5f\t' % float(self.ResultDict[result_property][time_t]))
            outfile.write('\n')


    def save_AdjMat(self, t):
        """
        :Purpose:
        Save Adjacency matrix in sparse format.
        """
        OutFileDir = os.path.expanduser(os.path.join(self.current_dir, 'Results'))
        if not os.path.isdir(OutFileDir):  # create directory if not existing
            os.mkdir(OutFileDir)
        OutFileName = os.path.join(OutFileDir, 'AdjacencyMatrix_atTime_%d.txt' % t)
        if os.path.isfile(OutFileName): os.remove(OutFileName)
        outfile = open(OutFileName, 'w')
        for n, row in enumerate(self.AdjMat.rows):
            outfile.write('%d:\t' % n)
            for partner in row:
                outfile.write('%d,' % partner)
            outfile.write('\n')
