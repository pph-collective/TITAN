#!/usr/bin/env python
# encoding: utf-8

from random import Random
import copy
from copy import deepcopy
import unittest
from scipy.stats import poisson

try:
    from .agent import *
except ImportError:
    raise ImportError("Can't import AgentClass")
from . import params


class PopulationClass:
    """
    :Purpose:
        This class constructs and represents the model population

    :Input:

        n : int
            Number of agents. Default: 10000

        r_seed : randomization seed

        model :str - one of "PrEP", "Incar", "NoIncar"

    :Attributes:

        :py:attr:`PopulationSize` : int
            Size of the population.

        :py:attr:`propIDU` : float
            Prevalence of intravenous drug users.

        :py:attr:`numIDU` : int
            Number of intravenous drug users.

        :py:attr:`propNIDU` : float
            Prevalence of non-intravenous drug users.

        :py:attr:`numNIDU` : int
            Number of non-intravenous drug users.

        :py:attr:`propND` : float
            Prevalence of not drug users.

        :py:attr:`numND` : int
            Number of not drug users.

        :py:attr:`Agents`: dict
            Dictionary of agents and their characteristics. The agents are
            the `keys' and a dicotionary of 'characteristic:value'
            pair is the entry.

        :py:attr:`IDU_agents`: list
            IDU drug users

        :py:attr:`NIDU_agents`: list
            NIDU drug users

        :py:attr:`ND_agents`: list
            ND drug users

        :py:attr:`MSM_agents`: list
            MSM agents

        :py:attr:`HIV_agents`: list
            HIV+ users

        :py:attr:`AIDS_agents`: list
            Users with AIDS

    :Methods:
        :py:meth:`get_agent_characteristic` \n
        :py:meth:`_set_population` \n
        :py:meth:`get_agents` \n
        :py:meth:`get_info_DrugUserType` \n
        :py:meth:`get_info_HIV_IDU` \n
        :py:meth:`get_info_DrugSexType`
    """

    def __init__(self, n=10000, rSeed=0, model=None):
        """
        :Purpose:
            Initialize PopulationClass object.
        """
        # Init RNG for population creation to rSeed
        self.popRandom = Random(rSeed)
        if type(n) is not int:
            raise ValueError("Population size must be integer")
        else:
            self.PopulationSize = n

        # Parameters
        self.numWhite = int(
            params.DemographicParams["WHITE"]["ALL"]["Proportion"] * self.PopulationSize
        )
        self.numBlack = int(
            params.DemographicParams["BLACK"]["ALL"]["Proportion"] * self.PopulationSize
        )
        # Nested dictionary for probability lookups by race
        # StratW = White, StratB = Black
        # HM incar @ 0.0279

        # MODEL TYPES
        MT_PrEP = False
        MT_Incar = False
        MT_NoIncar = False

        if model == "PrEP":
            MT_PrEP = True
        elif model == "Incar":
            MT_Incar = True
        elif model == "NoIncar":
            MT_NoIncar = True

        # Quick calc based on traj: 0.596 Female -> 0,01175/0.75 HIV -> 0.75 Tested+ = 11.75% HIV DIAG
        StratW = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}}
        if MT_PrEP:
            StratW["MSM"] = {
                "POP": 1.0,
                "HIV": 0.036195675,
                "AIDS": 0.51,
                "HAARTprev": 0.33,
                "INCARprev": 0.0,
                "TestedPrev": 0.816,
                "mNPart": 3,
            }
            StratW["HM"] = {
                "POP": 0.0,
                "HIV": 0.0,
                "AIDS": 0.0,
                "HAARTprev": 0.0,
                "INCARprev": 0.0,
                "TestedPrev": 0.0,
            }
            StratW["HF"] = {
                "POP": 0.0,
                "HIV": 0.0,
                "AIDS": 0.0,
                "HAARTprev": 0.0,
                "INCARprev": 0.0,
                "TestedPrev": 0.0,
            }
            StratW["IDU"] = {
                "POP": 0.0,
                "HIV": 0.0,
                "AIDS": 0.0,
                "HAARTprev": 0.0,
                "INCARprev": 0.0,
                "TestedPrev": 0.0,
            }
        elif MT_Incar:
            StratW["MSM"] = {
                "POP": 0.0,
                "HIV": 0.05,
                "AIDS": 0.05,
                "HAARTprev": 0.0,
                "INCARprev": 0.0,
                "TestedPrev": 0.86,
                "mNPart": 3,
            }
            StratW["HM"] = {
                "POP": 0.4044435412,
                "HIV": 0.0158 / 0.75,
                "AIDS": 0.678,
                "HAARTprev": 0.03748,
                "INCARprev": 0.0279,
                "TestedPrev": 0.75,
                "mNPart": 3,
            }
            StratW["HF"] = {
                "POP": 0.5955564588,
                "HIV": 0.01175 / 0.75,
                "AIDS": 0.573,
                "HAARTprev": 0.04054,
                "INCARprev": 0.0,
                "TestedPrev": 0.75,
                "mNPart": 2,
            }
            StratW["IDU"] = {
                "POP": 0.0,
                "HIV": 0.0,
                "AIDS": 0.0,
                "HAARTprev": 0.0,
                "INCARprev": 0.0,
                "TestedPrev": 0.0,
                "mNPart": 3,
            }
        elif MT_NoIncar:
            StratW["MSM"] = {
                "POP": 0.0,
                "HIV": 0.05,
                "AIDS": 0.05,
                "HAARTprev": 0.0,
                "INCARprev": 0.0,
                "TestedPrev": 0.86,
                "mNPart": 3,
            }
            StratW["HM"] = {
                "POP": 0.4044435412,
                "HIV": 0.0158 / 0.75,
                "AIDS": 0.678,
                "HAARTprev": 0.03748,
                "INCARprev": 0.0,
                "TestedPrev": 0.75,
                "mNPart": 3,
            }
            StratW["HF"] = {
                "POP": 0.5955564588,
                "HIV": 0.01175 / 0.75,
                "AIDS": 0.573,
                "HAARTprev": 0.04054,
                "INCARprev": 0.0,
                "TestedPrev": 0.75,
                "mNPart": 2,
            }
            StratW["IDU"] = {
                "POP": 0.0,
                "HIV": 0.0,
                "AIDS": 0.0,
                "HAARTprev": 0.0,
                "INCARprev": 0.0,
                "TestedPrev": 0.0,
                "mNPart": 3,
            }

        StratB = {"MSM": {}, "HM": {}, "HF": {}, "IDU": {}}
        StratB["MSM"] = {
            "POP": 0.0,
            "HIV": 0.0,
            "AIDS": 0.0,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.0,
        }
        StratB["HM"] = {
            "POP": 0.0,
            "HIV": 0.0,
            "AIDS": 0.0,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.0,
        }
        StratB["HF"] = {
            "POP": 0.0,
            "HIV": 0.0,
            "AIDS": 0.0,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.0,
        }
        StratB["IDU"] = {
            "POP": 0.0,
            "HIV": 0.0,
            "AIDS": 0.0,
            "HAARTprev": 0.0,
            "INCARprev": 0.0,
            "TestedPrev": 0.0,
        }

        self.ProbLookUp = {"WHITE": StratW, "BLACK": StratB}
        # drug user prevalence (proportion)

        self.propIDU = 0  ##190/10000.0
        self.numIDU = int(self.propIDU * self.PopulationSize)
        self.propNIDU = 0  ##640/10000.0
        self.numNIDU = int(self.propNIDU * self.PopulationSize)
        self.propND = 0.995  ##9170/10000.0
        self.numND = self.PopulationSize - self.numIDU - self.numNIDU
        self.Agents = dict()  # Main Dictionary Agents #REVIEW - delete, update anything if need be

        # List of IDU, NIDU, NDs
        # shuffle all Agents
        allAgents = list(range(self.PopulationSize))

        print("\tBuilding class sets")
        self.totalAgentClass = Agent_set(0, "TotalAgents") #REVIEW

        # All agent set list
        self.All_agentSet = Agent_set(0, "AllAgents")

        # HIV status agent sets
        self.HIV_agentSet = Agent_set(
            1, "HIV", parent=self.All_agentSet, numerator=self.All_agentSet
        )
        self.HIV_AIDS_agentSet = Agent_set(
            2, "AIDS", parent=self.HIV_agentSet, numerator=self.HIV_agentSet
        )

        # Drug use agent sets
        self.drugUse_agentSet = Agent_set(1, "DU", parent=self.All_agentSet)
        self.DU_NIDU_agentSet = Agent_set(2, "NIDU", parent=self.drugUse_agentSet)
        self.DU_IDU_agentSet = Agent_set(2, "IDU", parent=self.drugUse_agentSet)
        self.DU_NDU_agentSet = Agent_set(2, "NDU", parent=self.drugUse_agentSet)

        # Treatment agent sets
        self.treatment_agentSet = Agent_set(1, "Trtmt", parent=self.All_agentSet)
        self.Trt_Tstd_agentSet = Agent_set(
            2, "Testd", parent=self.treatment_agentSet, numerator=self.HIV_agentSet
        )
        self.Trt_PrEP_agentSet = Agent_set(2, "PrEP", parent=self.treatment_agentSet)
        self.Trt_PrEPelig_agentSet = Agent_set(2, "PrePelig", parent=self.treatment_agentSet)
        self.Trt_ART_agentSet = Agent_set(
            2, "ART", parent=self.treatment_agentSet, numerator=self.HIV_agentSet
        )
        self.Trt_SNE_agentSet = Agent_set(2, "SNE", parent=self.treatment_agentSet)
        self.Trt_MAT_agentSet = Agent_set(2, "MAT", parent=self.treatment_agentSet)

        # Sexual orientation agent sets
        self.SO_agentSet = Agent_set(1, "SO", parent=self.All_agentSet, numerator=self.All_agentSet)
        self.SO_HF_agentSet = Agent_set(
            2, "HF", parent=self.SO_agentSet, numerator=self.SO_agentSet
        )
        self.SO_HM_agentSet = Agent_set(
            2, "HM", parent=self.SO_agentSet, numerator=self.SO_agentSet
        )
        self.SO_MSM_agentSet = Agent_set(
            2, "MSM", parent=self.SO_agentSet, numerator=self.SO_agentSet
        )
        self.SO_MSW_agentSet = Agent_set(
            2, "MSW", parent=self.SO_agentSet, numerator=self.SO_agentSet
        )

        # Racial agent sets
        self.racial_agentSet = Agent_set(1, "Race", parent=self.All_agentSet)
        self.Race_WHITE_agentSet = Agent_set(2, "WHITE", parent=self.racial_agentSet)
        self.Race_BLACK_agentSet = Agent_set(2, "BLACK", parent=self.racial_agentSet)

        # Incarcerated agent sets
        self.incarcerated_agentSet = Agent_set(1, "Incar", parent=self.All_agentSet)

        # High risk agent sets
        self.highrisk_agentsSet = Agent_set(1, "HRisk", parent=self.All_agentSet)


        self.Relationships = Agent_set(0, "Relationships")

        # Nested dictionary for probability look-up
        IDU = {"MSM": {}, "HM": {}, "WSW": {}, "HF": {}}
        IDU["MSM"] = {"HIV": 0.55, "AIDS": 0.058, "HAARTa": 0}
        IDU["HM"] = {"HIV": 0.42, "AIDS": 0.058, "HAARTa": 0}
        IDU["WSW"] = {"HIV": 0.53, "AIDS": 0.058, "HAARTa": 0}
        IDU["HF"] = {"HIV": 0.39, "AIDS": 0.058, "HAARTa": 0}
        NIDU = {"MSM": {}, "HM": {}, "WSW": {}, "HF": {}}
        NIDU["MSM"] = {"HIV": 0.18, "AIDS": 0.02, "HAARTa": 0}
        NIDU["HM"] = {"HIV": 0.048, "AIDS": 0.002, "HAARTa": 0}
        NIDU["WSW"] = {"HIV": 0.048, "AIDS": 0.002, "HAARTa": 0}
        NIDU["HF"] = {"HIV": 0.048, "AIDS": 0.002, "HAARTa": 0}
        ND = {"MSM": {}, "HM": {}, "WSW": {}, "HF": {}}
        ND["Type"] = ([0, 1, 2, 3], [0.469, 0.493, 0.022, 0.016])
        ND["MSM"] = {"HIV": 0.08, "AIDS": 0.02, "HAARTa": 0}
        ND["HM"] = {"HIV": 0.015, "AIDS": 0.0003, "HAARTa": 0}
        ND["WSW"] = {"HIV": 0.012, "AIDS": 0.0003, "HAARTa": 0}
        ND["HF"] = {"HIV": 0.012, "AIDS": 0.0003, "HAARTa": 0}

        # Create agents in allAgents list
        self.White_agents = deepcopy(allAgents[0 : self.numWhite])
        self.Black_agents = deepcopy(allAgents[self.numWhite :])

        print("\tCreating agents")

        for agent in self.White_agents:
            agent_cl = self._return_new_Agent_class(agent, "WHITE")
            self.create_agent(agent_cl, "WHITE")
        for agent in self.Black_agents:
            agent_cl = self._return_new_Agent_class(agent, "BLACK")
            self.create_agent(agent_cl, "BLACK")
        # jail stock duration?
        jailDuration = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}}
        jailDuration[1] = {"p_value": (0.14), "min": 1, "max": 13}
        jailDuration[2] = {"p_value": (0.09), "min": 13, "max": 26}
        jailDuration[3] = {"p_value": (0.20), "min": 26, "max": 78}
        jailDuration[4] = {"p_value": (0.11), "min": 78, "max": 130}
        jailDuration[5] = {"p_value": (0.16), "min": 130, "max": 260}
        jailDuration[6] = {"p_value": (0.30), "min": 260, "max": 520}
        prob_Incarc = params.DemographicParams["WHITE"]["HM"]["mNPart"]
        for tmpA in self.All_agentSet._members:

            dice = self.popRandom.random()
            prob_Incarc = params.DemographicParams[tmpA._race][tmpA._SO]["INCARprev"]
            if dice < prob_Incarc:
                toss = 2 #REVIEW why is this hardcoded? - sarah to review
                if toss == 1:  # JAIL
                    tmpA._incar_bool = True
                    tmpA._incar_time = int(self.popRandom.triangular(1, 9, 3))
                else:  # PRISON
                    tmpA._incar_bool = True
                    durationBin = current_p_value = 0
                    p = self.popRandom.random()
                    while p > current_p_value:
                        durationBin += 1
                        current_p_value += jailDuration[durationBin]["p_value"]
                    timestay = self.popRandom.randint(
                        jailDuration[durationBin]["min"], jailDuration[durationBin]["max"]
                    )
                self.incarcerated_agentSet.add_agent(tmpA)


    def _return_agent_set(self):
        return self.totalAgentClass

    def _return_new_agent_dict(self, Deliminator, SexType="NULL"):
        """
        :Purpose:
        Return random agent dict of a new agent..
            Each agent is a key to an associated dictionary which stores the internal
            characteristics in form of an additinoal dictionary of the form
            ``characteristic:value``.
        :Input:
        DrugType : str
            Either 'IDU','NIDU' or 'ND'
        :Output:
             agent_dict : dict
        """
        Drugtype = "NULL"

        # Determine sextype
        tmp_rnd = self.popRandom.random()
        while SexType == "NULL":
            if tmp_rnd < params.DemographicParams[Deliminator]["HM"]["POP"]:
                SexType = "HM"
            elif tmp_rnd < (
                params.DemographicParams[Deliminator]["HM"]["POP"]
                + params.DemographicParams[Deliminator]["HF"]["POP"]
            ):
                SexType = "HF"
            else:
                SexType = "MSM"

        tmp_rnd = self.popRandom.random()

        # Determine drugtype
        tmp_rnd = self.popRandom.random()

        if tmp_rnd < params.DemographicParams[Deliminator]["IDU"]["POP"]:
            DrugType = "IDU"
        else:
            DrugType = "ND"

        # HIV
        if DrugType == "IDU":
            prob_HIV = params.DemographicParams[Deliminator]["IDU"]["HIV"]
        else:
            prob_HIV = params.DemographicParams[Deliminator][SexType]["HIV"]

        if self._MSMW:
            prob_HIV *= 2

        if self.popRandom.random() < prob_HIV:
            HIVStatus = 1

            # if HIV AIDS possible
            if DrugType == "IDU":
                prob_AIDS = params.DemographicParams[Deliminator]["IDU"]["AIDS"]
            else:
                prob_AIDS = params.DemographicParams[Deliminator][SexType]["AIDS"]

            if self.popRandom.random() < prob_AIDS:
                AIDSStatus = 1
            else:
                AIDSStatus = 0

            # HIV testing params
            if DrugType == "IDU":
                prob_Tested = params.DemographicParams[Deliminator]["IDU"]["TestedPrev"]
            else:
                prob_Tested = params.DemographicParams[Deliminator][SexType]["TestedPrev"]

            if self.popRandom.random() < prob_Tested:
                TestedStatus = 1

                # if tested HAART possible
                if DrugType == "IDU":
                    prob_HAART = params.DemographicParams[Deliminator]["IDU"]["HAARTprev"]
                else:
                    prob_HAART = params.DemographicParams[Deliminator][SexType]["HAARTprev"]

                if self.popRandom.random() < prob_HAART:
                    HAARTStatus = 1
                else:
                    HAARTStatus = 0

            else:
                TestedStatus = 0
                HAARTStatus = 0

            # if HIV, how long has the agent had it? Random sample
            HIV_time = self.popRandom.randint(1, 42)

        else:
            HIVStatus = 0
            AIDSStatus = 0
            HAARTStatus = 0
            HIV_time = 0
            TestedStatus = 0

            if params.PrEP_startT == 0:
                prob_PrEP = params.PrEP_Target
            else:
                prob_PrEP = 0.0

            if self.popRandom.random() < prob_PrEP:
                PrEP_Status = 1
            else:
                PrEP_Status = 0

        # Incarceration
        if DrugType == "IDU":
            prob_Incarc = params.DemographicParams[Deliminator]["IDU"]["INCARprev"]
        else:
            prob_Incarc = params.DemographicParams[Deliminator][SexType]["INCARprev"]

        if self.popRandom.random() < prob_Incarc:
            toss = self.popRandom.choice((1, 2))
            if toss == 1:  # JAIL
                incar_time = int(self.popRandom.triangular(6, 21, 15))
            else:  # PRISON
                incar_time = int(self.popRandom.triangular(30, 96, 60))
        else:
            incar_time = 0

        agent_dict = {
            "Race": Deliminator,
            "Drug Type": DrugType,
            "Sex Type": SexType,
            "HIV": HIVStatus,
            "Tested": TestedStatus,
            "AIDS": AIDSStatus,
            "HAARTa": HAARTStatus,
            "PrEP": PrEP_Status,
            "incar_t": incar_time,
            "HIV_t": HIV_time,
        }

        return agent_dict


    def _return_new_Agent_class(self, agentID, Race, SexType="NULL"):
        """
        :Purpose:
        Return random agent dict of a new agent..
            Each agent is a key to an associated dictionary which stores the internal
            characteristics in form of an additinoal dictionary of the form
            ``characteristic:value``.
        :Input:
        DrugType : str
            Either 'IDU','NIDU' or 'ND'
        :Output:
             agent_dict : dict
        """
        Drugtype = "NULL"

        # Determine sextype
        demBinP = 0.0
        tmp_rnd = self.popRandom.random()
        while SexType == "NULL":
            # For each demographic class within race
            for demClass in list(params.RaceClass1.keys()):
                # If demClass is enabled for model
                if demClass in params.agentSexTypes:
                    # Calculate probability
                    demBinP += params.DemographicParams[Race][demClass]["POP"]
                    # If match, set SexType
                    if tmp_rnd < demBinP:
                        SexType = demClass
                        break

        # Determine drugtype
        tmp_rnd = self.popRandom.random()

        # todo: FIX THIS TO GET BACK IDU
        if tmp_rnd < params.DemographicParams[Race]["IDU"]["POP"]:
            DrugType = "IDU"
        else:
            DrugType = "NDU"

        age, ageBin = self.getAge(Race)

        newAgent = Agent(agentID, SexType, age, Race, DrugType)
        newAgent._ageBin = ageBin

        if params.setting == 'Phil2005' and SexType == 'HM':
            if tmp_rnd < 0.06:
                newAgent._MSMW = True

        # HIV
        if DrugType == "IDU":
            prob_HIV = params.DemographicParams[Race]["IDU"]["HIV"]
        else:
            prob_HIV = params.DemographicParams[Race][SexType]["HIV"]

        if self.popRandom.random() < prob_HIV:
            HIVStatus = 1
            newAgent._HIV_bool = True

            # if HIV AIDS possible
            if DrugType == "IDU":
                prob_AIDS = params.DemographicParams[Race]["IDU"]["AIDS"]
            else:
                prob_AIDS = params.DemographicParams[Race][SexType]["AIDS"]

            if self.popRandom.random() < prob_AIDS:
                AIDSStatus = 1
                newAgent._AIDS_bool = True
            else:
                AIDSStatus = 0

            # HIV testing params
            if DrugType == "IDU":
                prob_Tested = params.DemographicParams[Race]["IDU"]["TestedPrev"]
            else:
                prob_Tested = params.DemographicParams[Race][SexType]["TestedPrev"]

            if self.popRandom.random() < prob_Tested:
                TestedStatus = 1
                newAgent._tested = True

                # if tested HAART possible
                if DrugType == "IDU":
                    prob_HAART = params.DemographicParams[Race]["IDU"]["HAARTprev"]
                else:
                    prob_HAART = params.DemographicParams[Race][SexType]["HAARTprev"]

                if self.popRandom.random() < prob_HAART:
                    HAARTStatus = 1
                    newAgent._HAART_bool = True
                    newAgent._treatment_bool = True
                else:
                    HAARTStatus = 0

            else:
                TestedStatus = 0
                HAARTStatus = 0

            # if HIV, how long has the agent had it? Random sample
            newAgent._HIV_time = self.popRandom.randint(1, 42)

        else:
            HIVStatus = 0
            AIDSStatus = 0
            HAARTStatus = 0
            HIV_time = 0
            TestedStatus = 0

            if params.flag_PrEP:
                if params.PrEP_startT == -1:
                    prob_PrEP = params.PrEP_Target
                else:
                    prob_PrEP = 0.0

                if self.popRandom.random() < prob_PrEP:
                    PrEP_Status = 1
                    newAgent._PrEP_bool = True
                    newAgent._treatment_bool = True
                else:
                    PrEP_Status = 0

        # Check if agent is HR as baseline.
        if self.popRandom.random() < params.DemographicParams[Race][SexType]["HighRiskPrev"]:
            newAgent._highrisk_bool = True
            newAgent._everhighrisk_bool = True

        diceroll = self.popRandom.random()

        if DrugType == 'IDU':
            if diceroll < 0.389:
                mNPart = 1
            elif diceroll < 0.389 + 0.150:
                mNPart = 2
            elif diceroll < 0.389 + 0.150 + 0.089:
                mNPart = 3
            elif diceroll < 0.389 + 0.150 + 0.089 + 0.067:
                mNPart = 4
            elif diceroll < 0.389 + 0.150 + 0.089 + 0.067 + 0.098:
                mNPart = self.popRandom.randrange(5, 6, 1)
            elif diceroll < 0.389 + 0.150 + 0.089 + 0.067 + 0.098 + 0.108:
                mNPart = self.popRandom.randrange(7, 10, 1)
            elif diceroll < 0.389 + 0.150 + 0.089 + 0.067 + 0.098 + 0.108 + 0.02212:
                mNPart = self.popRandom.randrange(17, 30, 1)
            else:
                mNPart = self.popRandom.randrange(31, 60, 1)
        else:
            if diceroll < 0.84:
                mNPart = 1
            elif diceroll < 0.84 + 0.13:
                mNPart = 2
            else:
                mNPart = self.popRandom.randrange(3, 4, 1)

        # Partnership demographics
        if params.setting == "Scott":
            newAgent._mean_num_partners = (
                mNPart
            )
        else:
            newAgent._mean_num_partners = poisson.rvs(
                params.DemographicParams[Race][SexType]["NUMPartn"], size=1
            )

        return newAgent


    def create_agent(self, agent_cl, Deliminator):
        """
        :Purpose:
            Creat a new agent in the population.
                Each agent is a key to an associated dictionary which stores the internal
                characteristics in form of an additinoal dictionary of the form
                ``characteristic:value``.

        :Input:
            agent : int

            Deliminator : str currently race

        """

        def addToSubsets(targetSet, agent, agentParam=None):
            try:
                targetSet.add_agent(agent)
            except:
                print(
                    "agent %s is already a member of agent set %s"
                    % (agent.get_ID(), targetSet.get_ID())
                )

            if agentParam:
                try:
                    targetSet._subset[agentParam].add_agent(agent)
                except:
                    print(
                        "agent %s is already a member of agent set %s"
                        % (agent_cl.get_ID(), targetSet._subset[agentParam].get_ID())
                    )

        # Add to all agent set
        self.All_agentSet.add_agent(agent_cl)

        # Add to correct SO set
        addToSubsets(self.SO_agentSet, agent_cl, agent_cl._SO)

        # Add to correct DU set
        addToSubsets(self.drugUse_agentSet, agent_cl, agent_cl._DU)

        # Add to correct racial set
        addToSubsets(self.racial_agentSet, agent_cl, agent_cl._race)

        if agent_cl._HIV_bool:
            addToSubsets(self.HIV_agentSet, agent_cl)
            if agent_cl._AIDS_bool:
                addToSubsets(self.HIV_AIDS_agentSet, agent_cl)

        # Add to correct treatment set
        if agent_cl._treatment_bool:
            addToSubsets(self.treatment_agentSet, agent_cl)
            if agent_cl._HAART_bool:
                addToSubsets(self.Trt_ART_agentSet, agent_cl)

        if agent_cl._PrEP_bool:
            addToSubsets(self.Trt_PrEP_agentSet, agent_cl)

        if agent_cl._tested:
            addToSubsets(self.Trt_Tstd_agentSet, agent_cl)

        if agent_cl._incar_bool:
            addToSubsets(self.incarcerated_agentSet, agent_cl)

        if agent_cl._highrisk_bool:
            addToSubsets(self.highrisk_agentsSet, agent_cl)


    def getAge(self, race):
        rand = self.popRandom.random()
        minAge = 15
        maxAge = 80
        ageBin = 0

        if params.setting == "AtlantaMSM":
            if rand < params.ageMatrix[race]["Prop"][1]:
                minAge = 18
                maxAge = 19
                ageBin = 1
            elif rand < params.ageMatrix[race]["Prop"][2]:
                minAge = 20
                maxAge = 24
                ageBin = 2
            elif rand < params.ageMatrix[race]["Prop"][3]:
                minAge = 25
                maxAge = 29
                ageBin = 3
            elif rand <= params.ageMatrix[race]["Prop"][4]:
                minAge = 30
                maxAge = 39
                ageBin = 4
        else:
            if rand < params.ageMatrix[race]["Prop"][1]:
                minAge = 15
                maxAge = 24
                ageBin = 1
            elif rand < params.ageMatrix[race]["Prop"][2]:
                minAge = 25
                maxAge = 34
                ageBin = 2
            elif rand < params.ageMatrix[race]["Prop"][3]:
                minAge = 35
                maxAge = 44
                ageBin = 3
            elif rand < params.ageMatrix[race]["Prop"][4]:
                minAge = 45
                maxAge = 54
                ageBin = 4
            else:
                minAge = 55
                maxAge = 80
                ageBin = 5

        age = self.popRandom.randrange(minAge, maxAge)
        return age, ageBin
