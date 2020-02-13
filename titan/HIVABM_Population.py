#!/usr/bin/env python
# encoding: utf-8

from random import Random
from copy import deepcopy
from typing import Sequence, List, Dict, Optional, Any
from scipy.stats import poisson  # type: ignore
import numpy as np  # type: ignore
from dotmap import DotMap  # type: ignore

from .agent import Agent_set, Agent, Relationship
from .ABM_partnering import get_partner, get_partnership_duration
from . import probabilities as prob


class PopulationClass:
    """
    :Purpose:
        This class constructs and represents the model population

    :Input:

        n : int
            Number of agents

        r_seed : randomization seed


        params : Map containing parameters

    """

    def __init__(self, pop_seed: int, params: DotMap):
        """
        :Purpose:
            Initialize PopulationClass object.
        """
        # Init RNG for population creation to pop_seed
        self.pop_random = Random(pop_seed)
        np.random.seed(pop_seed)

        self.params = params

        # build weights of population sex types by race - SARAH READ THIS
        self.pop_weights: Dict[str, Dict[str, List[Any]]] = {}
        for race in params.classes.races:
            self.pop_weights[race] = {}
            self.pop_weights[race]["values"] = []
            self.pop_weights[race]["weights"] = []
            for st in params.classes.sex_types:
                self.pop_weights[race]["values"].append(st)
                self.pop_weights[race]["weights"].append(
                    params.demographics[race][st].ppl
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
        self.treatment_agentSet = Agent_set("Trtmt", parent=self.All_agentSet)
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
        self.SO_MTF_agentSet = Agent_set(
            "MTF", parent=self.SO_agentSet, numerator=self.SO_agentSet
        )

        self.SO_WSW_agentSet = Agent_set(
            "WSW", parent=self.SO_agentSet, numerator=self.SO_agentSet
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

        for race in params.classes.races:
            for i in range(round(params.model.num_pop * params.demographics[race].ppl)):
                agent = self.create_agent(race)
                self.add_agent_to_pop(agent)

        self.initialize_incarceration()

    def initialize_incarceration(self):
        jail_duration = prob.jail_duration()

        for a in self.All_agentSet._members:

            prob_incar = self.params.demographics[a._race][a._SO].incar.init
            if self.pop_random.random() < prob_incar:
                a._incar_bool = True
                bin = current_p_value = 0
                p = self.pop_random.random()

                while p > current_p_value:
                    bin += 1
                    current_p_value += jail_duration[bin]["p_value"]

                self.incarcerated_agentSet.add_agent(a)

    def create_agent(self, race: str, sex_type: str = "NULL") -> Agent:
        """
        :Purpose:
            Return a new agent according to population characteristics
        :Input:
            race : "BLACK" or "WHITE"
            sex_type : default "NULL"
        :Output:
             agent : Agent
        """
        if sex_type == "NULL":
            sex_type = np.random.choice(
                self.pop_weights[race]["values"], p=self.pop_weights[race]["weights"]
            )

        # Determine drugtype
        # todo: FIX THIS TO GET BACK IDU
        if self.pop_random.random() < self.params.demographics[race]["IDU"].ppl:
            drug_type = "IDU"
        else:
            drug_type = "NDU"

        age, age_bin = self.get_age(race)

        agent = Agent(sex_type, age, race, drug_type)
        agent._ageBin = age_bin

        # TO_REVIEW setting based logic
        # if params.setting == "Phil2005" and sex_type == "HM":
        #     if self.pop_random.random() < 0.06:
        #         agent._MSMW = True

        if drug_type == "IDU":
            agent_params = self.params.demographics[race]["IDU"]
        else:
            agent_params = self.params.demographics[race][sex_type]

        # HIV
        if self.pop_random.random() < agent_params.hiv.init:
            agent._HIV_bool = True

            if self.pop_random.random() < agent_params.aids.init:
                agent._AIDS_bool = True

            if self.pop_random.random() < agent_params.hiv.dx.init:
                agent._tested = True

                if self.pop_random.random() < agent_params.haart.init:
                    agent._HAART_bool = True
                    agent._treatment_bool = True

                    haart_adh = self.params.demographics[race][sex_type].haart.adherence
                    if self.pop_random.random() < haart_adh:
                        adherence = 5
                    else:
                        adherence = self.pop_random.randint(1, 4)

                    # add to agent haart set
                    agent._HAART_adh = adherence
                    agent._HAART_time = 0

            # if HIV, how long has the agent had it? Random sample
            agent._HIV_time = self.pop_random.randint(1, 42)

        else:

            if self.params.features.prep:
                if self.params.prep.start == 0:
                    prob_prep = self.params.prep.target
                else:
                    prob_prep = 0.0

                if self.pop_random.random() < prob_prep:
                    agent._PrEP_bool = True
                    agent._treatment_bool = True

        # Check if agent is HR as baseline.
        if (
            self.pop_random.random()
            < self.params.demographics[race][sex_type].high_risk.init
        ):
            agent._highrisk_bool = True
            agent._everhighrisk_bool = True

        # Partnership demographics # TO_REVIEW setting specific logic
        # if params.setting == "Scott":
        #     agent._mean_num_partners = prob.get_mean_num_partners(
        #         drug_type, self.pop_random
        #     )
        # else:
        agent._mean_num_partners = poisson.rvs(
            self.params.demographics[race][sex_type].num_partners, size=1
        )

        return agent

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

        def add_to_subsets(target, agent, agent_param=None):
            target.add_agent(agent)
            if agent_param:
                target._subset[agent_param].add_agent(agent)

        # Add to all agent set
        self.All_agentSet.add_agent(agent)

        # Add to correct SO set
        add_to_subsets(self.SO_agentSet, agent, agent._SO)

        # Add to correct DU set
        add_to_subsets(self.drugUse_agentSet, agent, agent._DU)

        # Add to correct racial set
        add_to_subsets(self.racial_agentSet, agent, agent._race)

        if agent._HIV_bool:
            add_to_subsets(self.HIV_agentSet, agent)
            if agent._AIDS_bool:
                add_to_subsets(self.HIV_AIDS_agentSet, agent)

        # Add to correct treatment set
        if agent._treatment_bool:
            add_to_subsets(self.treatment_agentSet, agent)
            if agent._HAART_bool:
                add_to_subsets(self.Trt_ART_agentSet, agent)

        if agent._PrEP_bool:
            add_to_subsets(self.Trt_PrEP_agentSet, agent)

        if agent._tested:
            add_to_subsets(self.Trt_Tstd_agentSet, agent)

        if agent._incar_bool:
            add_to_subsets(self.incarcerated_agentSet, agent)

        if agent._highrisk_bool:
            add_to_subsets(self.highrisk_agentsSet, agent)

    def get_age(self, race: str):
        rand = self.pop_random.random()

        # REVIEWED why does AtlantaMSM use different age bins? should this all be paramable? - this will be revisited with future age things - got rid of it here, but Atlanta's setting's demograhpics.race.rage needs to be updated to match what was here

        bins = self.params.demographics[race].age

        for i in range(1, 6):
            if rand < bins[i].prob:
                min_age = bins[i].min
                max_age = bins[i].max
                break

        age = self.pop_random.randrange(min_age, max_age)
        return age, i

    # REVIEWED should these be in the network class? - max to incorporate with network/pop/model disentangling?
    def update_agent_partners(self, graph, agent: Agent, params: DotMap):
        partner = get_partner(agent, self.All_agentSet, params)

        if partner:
            duration = get_partnership_duration(agent, params)
            rel = Relationship(agent, partner, duration)
            self.Relationships.append(rel)
            graph.add_edge(rel._ID1, rel._ID2)
        else:
            graph.add_node(agent)

    def update_partner_assignments(self, graph, params: DotMap):
        # Now create partnerships until available partnerships are out
        eligible_agents = self.All_agentSet
        for agent in eligible_agents.iter_agents():
            acquirePartnerProb = params.calibration.sex.partner * (
                agent._mean_num_partners / (12.0)
            )
            if self.pop_random.random() < acquirePartnerProb:
                self.update_agent_partners(graph, agent, params)
            else:
                graph.add_node(
                    agent
                )  # TO_REVIEW shouldn't adding the agent as a node happen no matter what?
