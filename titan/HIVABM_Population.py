#!/usr/bin/env python
# encoding: utf-8

from random import Random
from copy import deepcopy
from typing import Sequence, List, Dict, Optional, Any
from scipy.stats import poisson  # type: ignore
import numpy as np  # type: ignore

from .agent import Agent_set, Agent, Relationship
from .ABM_partnering import get_partner, get_partnership_duration
from . import params  # type: ignore
from . import probabilities as prob


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

    def __init__(self, n: int = 10000, rSeed: int = 0, model: str = None):
        """
        :Purpose:
            Initialize PopulationClass object.
        """
        # Init RNG for population creation to rSeed
        self.popRandom = Random(rSeed)
        np.random.seed(rSeed)

        if type(n) is not int:
            raise ValueError("Population size must be integer")
        else:
            self.PopulationSize = (
                n
            )  # REVIWED PopulationSize not really used and is calculable - needed? - nope

        # Parameters
        self.numWhite = round(
            params.DemographicParams["WHITE"]["ALL"]["Proportion"] * self.PopulationSize
        )
        self.numBlack = round(
            params.DemographicParams["BLACK"]["ALL"]["Proportion"] * self.PopulationSize
        )

        # build weights of population sex types by race - SARAH READ THIS
        self.pop_weights: Dict[str, Dict[str, List[Any]]] = {}
        for race in params.DemographicParams.keys():
            self.pop_weights[race] = {}
            self.pop_weights[race]["values"] = []
            self.pop_weights[race]["weights"] = []
            for st in params.agentSexTypes:
                self.pop_weights[race]["values"].append(st)
                self.pop_weights[race]["weights"].append(
                    params.DemographicParams[race][st]["POP"]
                )

        print("\tBuilding class sets")

        # All agent set list
        self.All_agentSet = Agent_set("AllAgents")

        # HIV status agent sets
        self.HIV_agentSet = Agent_set(
            "HIV", parent=self.All_agentSet, numerator=self.All_agentSet
        )
        self.HIV_AIDS_agentSet = Agent_set(
            "AIDS", parent=self.HIV_agentSet, numerator=self.HIV_agentSet
        )

        # Drug use agent sets
        self.drugUse_agentSet = Agent_set("DU", parent=self.All_agentSet)
        self.DU_NIDU_agentSet = Agent_set("NIDU", parent=self.drugUse_agentSet)
        self.DU_IDU_agentSet = Agent_set("IDU", parent=self.drugUse_agentSet)
        self.DU_NDU_agentSet = Agent_set("NDU", parent=self.drugUse_agentSet)

        # Treatment agent sets
        self.aware_agentSet = Agent_set("LAI aware", parent=self.All_agentSet)
        self.treatment_agentSet = Agent_set("Trtmt", parent=self.All_agentSet)
        self.PCA_agentSet = Agent_set("PCA", parent=self.All_agentSet)
        self.Trt_Tstd_agentSet = Agent_set(
            "Testd", parent=self.treatment_agentSet, numerator=self.HIV_agentSet
        )
        self.Trt_PrEP_agentSet = Agent_set("PrEP", parent=self.treatment_agentSet)
        self.Trt_PrEPelig_agentSet = Agent_set(
            "PrePelig", parent=self.treatment_agentSet
        )
        self.Trt_ART_agentSet = Agent_set(
            "ART", parent=self.treatment_agentSet, numerator=self.HIV_agentSet
        )

        # Sexual orientation agent sets
        self.SO_agentSet = Agent_set(
            "SO", parent=self.All_agentSet, numerator=self.All_agentSet
        )
        self.SO_HF_agentSet = Agent_set(
            "HF", parent=self.SO_agentSet, numerator=self.SO_agentSet
        )
        self.SO_HM_agentSet = Agent_set(
            "HM", parent=self.SO_agentSet, numerator=self.SO_agentSet
        )
        self.SO_MSM_agentSet = Agent_set(
            "MSM", parent=self.SO_agentSet, numerator=self.SO_agentSet
        )
        self.SO_MSW_agentSet = Agent_set(
            "MSW", parent=self.SO_agentSet, numerator=self.SO_agentSet
        )

        # Racial agent sets
        self.racial_agentSet = Agent_set("Race", parent=self.All_agentSet)
        self.Race_WHITE_agentSet = Agent_set("WHITE", parent=self.racial_agentSet)
        self.Race_BLACK_agentSet = Agent_set("BLACK", parent=self.racial_agentSet)

        # Incarcerated agent sets
        self.incarcerated_agentSet = Agent_set("Incar", parent=self.All_agentSet)

        # High risk agent sets
        self.highrisk_agentsSet = Agent_set("HRisk", parent=self.All_agentSet)

        self.Relationships: List[Relationship] = []

        print("\tCreating agents")

        for i in range(self.numWhite):
            agent = self.create_agent("WHITE")
            self.add_agent_to_pop(agent)

        for i in range(self.numBlack):
            agent = self.create_agent("BLACK")
            self.add_agent_to_pop(agent)

        self.initialize_incarceration()

    def initialize_incarceration(self):
        jailDuration = prob.jail_duration()

        for tmpA in self.All_agentSet._members:

            prob_Incarc = params.DemographicParams[tmpA._race][tmpA._SO]["INCARprev"]
            if self.popRandom.random() < prob_Incarc:
                tmpA._incar_bool = True
                durationBin = current_p_value = 0
                p = self.popRandom.random()

                while p > current_p_value:
                    durationBin += 1
                    current_p_value += jailDuration[durationBin]["p_value"]

                self.incarcerated_agentSet.add_agent(tmpA)

    def create_agent(self, Race: str, SexType: str = "NULL") -> Agent:
        """
        :Purpose:
            Return a new agent according to population characteristics
        :Input:
            Race : "BLACK" or "WHITE"
            SexType : default "NULL"
        :Output:
             newAgent : Agent
        """
        if SexType == "NULL":
            SexType = np.random.choice(
                self.pop_weights[Race]["values"], p=self.pop_weights[Race]["weights"]
            )

        # Determine drugtype
        # todo: FIX THIS TO GET BACK IDU
        if self.popRandom.random() < params.DemographicParams[Race]["IDU"]["POP"]:
            DrugType = "IDU"
        elif self.popRandom.random() < params.DemographicParams[Race][SexType]["nidu"]:
            DrugType = "NIDU"
        else:
            DrugType = "NDU"

        age, ageBin = self.get_age(Race)

        newAgent = Agent(SexType, age, Race, DrugType)
        newAgent._ageBin = ageBin

        if params.setting == "Phil2005" and SexType == "HM":  # TODO: make this a flag
            if self.popRandom.random() < 0.06:
                newAgent._MSMW = True

        # HIV
        if DrugType == "IDU":
            prob_HIV = params.DemographicParams[Race]["IDU"]["HIV"]
        else:
            prob_HIV = params.DemographicParams[Race][SexType]["HIV"]

        if self.popRandom.random() < prob_HIV:
            newAgent._HIV_bool = True

            # if HIV AIDS possible
            if DrugType == "IDU":
                prob_AIDS = params.DemographicParams[Race]["IDU"]["AIDS"]
            else:
                prob_AIDS = params.DemographicParams[Race][SexType]["AIDS"]

            if self.popRandom.random() < prob_AIDS:
                newAgent._AIDS_bool = True

            # HIV testing params
            if DrugType == "IDU":
                prob_Tested = params.DemographicParams[Race]["IDU"]["TestedPrev"]
            else:
                prob_Tested = params.DemographicParams[Race][SexType]["TestedPrev"]

            if self.popRandom.random() < prob_Tested:
                newAgent._tested = True

                # if tested HAART possible
                if DrugType == "IDU":
                    prob_HAART = params.DemographicParams[Race]["IDU"]["HAARTprev"]
                else:
                    prob_HAART = params.DemographicParams[Race][SexType]["HAARTprev"]

                if self.popRandom.random() < prob_HAART:
                    newAgent._HAART_bool = True
                    newAgent._treatment_bool = True

                    HAART_ADH = params.DemographicParams[Race][SexType]["HAARTadh"]
                    if self.popRandom.random() < HAART_ADH:
                        adherence = 5
                    else:
                        adherence = self.popRandom.randint(1, 4)

                    # add to agent haart set
                    newAgent._HAART_adh = adherence
                    newAgent._HAART_time = 0

            # if HIV, how long has the agent had it? Random sample
            newAgent._HIV_time = self.popRandom.randint(1, 42)

        else:
            if params.flag_PrEP:
                if params.PrEP_startT == -1:
                    prob_PrEP = params.PrEP_Target
                else:
                    prob_PrEP = 0.0

                if self.popRandom.random() < prob_PrEP:
                    newAgent._PrEP_bool = True
                    newAgent._treatment_bool = True
                    if (
                        self.popRandom.random() > params.LAI_chance
                        and "Inj" in params.PrEP_type
                    ):
                        newAgent.PrEP_type = "Inj"
                    else:
                        newAgent.PrEP_type = "Oral"

        # Check if agent is HR as baseline.
        if (
            self.popRandom.random()
            < params.DemographicParams[Race][SexType]["HighRiskPrev"]
        ):
            newAgent._highrisk_bool = True
            newAgent._everhighrisk_bool = True

        # Partnership demographics
        if params.setting == "Scott":
            newAgent._mean_num_partners = prob.get_mean_num_partners(
                DrugType, self.popRandom
            )
        elif params.mean_partner_type == "bins":
            assert (
                sum(params.partnerNumber.values()) >= 1
            ), "Bin probabilities must add up to 1!"
            pn_prob = self.popRandom.random()
            current_p_value = ptnBin = 0

            while pn_prob > current_p_value:
                current_p_value += params.partnerNumber[ptnBin]
                ptnBin += 1
            newAgent._mean_num_partners = ptnBin
        else:
            newAgent._mean_num_partners = poisson.rvs(
                params.DemographicParams[Race][SexType]["NUMPartn"], size=1
            )

        if params.flag_PCA:
            if self.popRandom.random() < params.starting_awareness:
                newAgent.awareness = True
            attprob = self.popRandom.random()
            pvalue = 0.0
            for k, v in params.attitude.items():
                pvalue += v
                if attprob < pvalue:
                    newAgent.opinion = v
                    break
            assert newAgent.opinion in range(
                5
            ), (
                "Agents opinion of injectible PrEP is out of bounds"
            )  # TODO: move to testing

        return newAgent

    def add_agent_to_pop(self, agent: Agent):
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
            targetSet.add_agent(agent)
            if agentParam:
                targetSet._subset[agentParam].add_agent(agent)

        # Add to all agent set
        self.All_agentSet.add_agent(agent)

        # Add to correct SO set
        addToSubsets(self.SO_agentSet, agent, agent._SO)

        # Add to correct DU set
        addToSubsets(self.drugUse_agentSet, agent, agent._DU)

        # Add to correct racial set
        addToSubsets(self.racial_agentSet, agent, agent._race)

        if agent._HIV_bool:
            addToSubsets(self.HIV_agentSet, agent)
            if agent._AIDS_bool:
                addToSubsets(self.HIV_AIDS_agentSet, agent)

        # Add to correct treatment set

        if agent.awareness:
            addToSubsets(self.aware_agentSet, agent)

        if agent._PrEP_bool:
            addToSubsets(self.Trt_PrEP_agentSet, agent)
        if agent._treatment_bool:
            addToSubsets(self.treatment_agentSet, agent)
            if agent._HAART_bool:
                addToSubsets(self.Trt_ART_agentSet, agent)

        if agent._PrEP_bool:
            addToSubsets(self.Trt_PrEP_agentSet, agent)

        if agent._tested:
            addToSubsets(self.Trt_Tstd_agentSet, agent)

        if agent._incar_bool:
            addToSubsets(self.incarcerated_agentSet, agent)

        if agent._highrisk_bool:
            addToSubsets(self.highrisk_agentsSet, agent)

    def get_age(self, race: str):
        rand = self.popRandom.random()

        # REVIEWED why does AtlantaMSM use different age bins? should this all be paramable? - this will be revisited with future age things
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
                minAge = 15
                maxAge = 80
                ageBin = 0
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

    # REVIEWED should these be in the network class? - max to incorporate with network/pop/model disentangling?
    def update_agent_partners(self, graph, agent: Agent) -> bool:
        partner = get_partner(agent, self.All_agentSet)
        noMatch = False
        rel_type = ""
        bond_type = "sexOnly"

        def bondtype(bond_dict):
            pvalue = 0.0
            bond_probability = self.popRandom.random()
            bonded_type = "sexualOnly"
            for reltype, p in bond_dict.items():
                pvalue += p
                if bond_probability < pvalue:
                    bonded_type = reltype
                    break
            return bonded_type

        if partner:
            duration = get_partnership_duration(agent)
            if params.bond_type:
                if agent._DU == "IDU" and partner._DU == "IDU":
                    bond_type = bondtype(params.bond_type_probs_IDU)
                else:
                    bond_type = bondtype(params.bond_type_probs)

            tmp_relationship = Relationship(
                agent, partner, duration, rel_type=bond_type
            )

            self.Relationships.append(tmp_relationship)
            graph.add_edge(
                tmp_relationship._ID1, tmp_relationship._ID2, relationship=rel_type
            )
        else:
            graph.add_node(agent)
            noMatch = True
        return noMatch

    def update_partner_assignments(self, partnerTurnover: float, graph):
        # Now create partnerships until available partnerships are out
        EligibleAgents = self.All_agentSet
        for agent in EligibleAgents.iter_agents():
            acquirePartnerProb = (
                params.cal_SexualPartScaling
                * partnerTurnover
                * (agent._mean_num_partners / (12.0))
            )
            if self.popRandom.random() < acquirePartnerProb:
                self.update_agent_partners(graph, agent)
            else:
                graph.add_node(agent)
