#!/usr/bin/env python
# encoding: utf-8

import random

from typing import List, Dict, Any
from scipy.stats import poisson  # type: ignore
import numpy as np  # type: ignore
from dotmap import DotMap  # type: ignore

from .agent import AgentSet, Agent, Relationship
from .partnering import get_partner, get_partnership_duration


class Population:
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
            Initialize Population object.
        """
        # Init RNG for population creation to pop_seed
        self.pop_random = random.Random(pop_seed)
        random.seed(
            pop_seed
        )  # this sets the global random seed for the population generation phase, during model init it gets reset at the very end
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
        self.all_agents = AgentSet("AllAgents")

        # HIV status agent sets
        self.hiv_agents = AgentSet(
            "HIV", parent=self.all_agents, numerator=self.all_agents
        )
        self.hiv_aids_agents = AgentSet(
            "AIDS", parent=self.hiv_agents, numerator=self.hiv_agents
        )

        # Drug use agent sets
        self.drug_use_agents = AgentSet("DU", parent=self.all_agents)
        self.drug_use_noninj_agents = AgentSet(
            "NonInj", parent=self.drug_use_agents
        )  # TO_REVIEW not really used
        self.drug_use_inj_agents = AgentSet("Inj", parent=self.drug_use_agents)
        self.drug_use_none_agents = AgentSet(
            "None", parent=self.drug_use_agents
        )  # TO_REVIEW not really used

        # Treatment agent sets
        self.intervention_agents = AgentSet("Trtmt", parent=self.all_agents)
        self.intervention_dx_agents = AgentSet(
            "Testd", parent=self.intervention_agents, numerator=self.hiv_agents
        )
        self.intervention_prep_agents = AgentSet(
            "PrEP", parent=self.intervention_agents
        )
        self.intervention_prep_eligible_agents = AgentSet(
            "PrePelig", parent=self.intervention_agents
        )  # TO_REVIEW not used anywhere
        self.intervention_haart_agents = AgentSet(
            "ART", parent=self.intervention_agents, numerator=self.hiv_agents
        )

        # Sexual orientation agent sets
        self.sex_type_agents = AgentSet(
            "SO", parent=self.all_agents, numerator=self.all_agents
        )
        self.sex_type_HF_agents = AgentSet(
            "HF", parent=self.sex_type_agents, numerator=self.sex_type_agents
        )
        self.sex_type_HM_agents = AgentSet(
            "HM", parent=self.sex_type_agents, numerator=self.sex_type_agents
        )
        self.sex_type_MSM_agents = AgentSet(
            "MSM", parent=self.sex_type_agents, numerator=self.sex_type_agents
        )
        self.sex_type_MTF_agents = AgentSet(
            "MTF", parent=self.sex_type_agents, numerator=self.sex_type_agents
        )
        self.sex_type_WSW_agents = AgentSet(
            "WSW", parent=self.sex_type_agents, numerator=self.sex_type_agents
        )

        # Racial agent sets
        self.race_agents = AgentSet("Race", parent=self.all_agents)
        self.race_white_agents = AgentSet("WHITE", parent=self.race_agents)
        self.race_black_agents = AgentSet("BLACK", parent=self.race_agents)

        # Incarcerated agent sets
        self.incarcerated_agents = AgentSet("Incar", parent=self.all_agents)

        # High risk agent sets
        self.high_risk_agents = AgentSet("HRisk", parent=self.all_agents)

        self.relationships: List[Relationship] = []

        print("\tCreating agents")

        for race in params.classes.races:
            for i in range(round(params.model.num_pop * params.demographics[race].ppl)):
                agent = self.create_agent(race)
                self.add_agent_to_pop(agent)

        self.initialize_incarceration()

    def initialize_incarceration(self):

        for a in self.all_agents.members:
            jail_duration = self.params.demographics[a.race][a.so].incar.duration.init

            prob_incar = self.params.demographics[a.race][a.so].incar.init
            if self.pop_random.random() < prob_incar:
                a.incar = True
                bin = current_p_value = 0
                p = self.pop_random.random()

                while p > current_p_value:
                    bin += 1
                    current_p_value += jail_duration[bin].prob

                a.incar_time = self.pop_random.randrange(
                    jail_duration[bin].min, jail_duration[bin].max
                )
                self.incarcerated_agents.add_agent(a)

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
        # todo: FIX THIS TO GET BACK PWID
        if self.pop_random.random() < self.params.demographics[race]["PWID"].ppl:
            drug_type = "Inj"
        else:
            drug_type = "None"

        age, age_bin = self.get_age(race)

        agent = Agent(sex_type, age, race, drug_type)
        agent.age_bin = age_bin

        if self.params.features.msmw and sex_type == "HM":
            if self.pop_random.random() < 0.06:
                agent.msmw = True

        if drug_type == "Inj":
            agent_params = self.params.demographics[race]["PWID"]
        else:
            agent_params = self.params.demographics[race][sex_type]

        # HIV
        if self.pop_random.random() < agent_params.hiv.init:
            agent.hiv = True

            if self.pop_random.random() < agent_params.aids.init:
                agent.aids = True

            if self.pop_random.random() < agent_params.hiv.dx.init:
                agent.hiv_dx = True

                if self.pop_random.random() < agent_params.haart.init:
                    agent.haart = True
                    agent.intervention_ever = True

                    haart_adh = self.params.demographics[race][sex_type].haart.adherence
                    if self.pop_random.random() < haart_adh:
                        adherence = 5
                    else:
                        adherence = self.pop_random.randint(1, 4)

                    # add to agent haart set
                    agent.haart_adherence = adherence
                    agent.haart_time = 0

            # if HIV, how long has the agent had it? Random sample
            agent.hiv_time = self.pop_random.randint(1, 42)

        else:

            if self.params.features.prep:
                if self.params.prep.start == 0:
                    prob_prep = self.params.prep.target
                else:
                    prob_prep = 0.0

                if self.pop_random.random() < prob_prep:
                    agent.prep = True
                    agent.intervention_ever = True
                    if (
                        self.pop_random.random() > self.params.prep.lai.prob
                        and "Inj" in self.params.prep.type
                    ):
                        agent.prep_type = "Inj"
                    else:
                        agent.prep_type = "Oral"

        # Check if agent is HR as baseline.
        if (
            self.params.features.high_risk
            and self.pop_random.random()
            < self.params.demographics[race][sex_type].high_risk.init
        ):
            agent.high_risk = True
            agent.high_risk_ever = True

        # Partnership demographics
        if self.params.model.population.num_partners.type == "bins":
            pn_prob = self.pop_random.random()
            current_p_value = bin = 0

            while pn_prob > current_p_value:
                current_p_value += self.params.model.population.num_partners.bins[
                    bin
                ].prob
                bin += 1
            agent.mean_num_partners = bin
        else:
            agent.mean_num_partners = poisson.rvs(
                self.params.demographics[race][sex_type].num_partners, size=1
            )

        if self.params.features.pca:
            if self.pop_random.random() < self.params.prep.pca.awareness.init:
                agent.awareness = True
            attprob = self.pop_random.random()
            pvalue = 0.0
            for bin, fields in self.params.pca.attitude.items():
                pvalue += fields.prob
                if attprob < pvalue:
                    agent.opinion = bin
                    break

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
                target.subset[agent_param].add_agent(agent)

        # Add to all agent set
        self.all_agents.add_agent(agent)

        # Add to correct SO set
        add_to_subsets(self.sex_type_agents, agent, agent.so)

        # Add to correct DU set
        add_to_subsets(self.drug_use_agents, agent, agent.drug_use)

        # Add to correct racial set
        add_to_subsets(self.race_agents, agent, agent.race)

        if agent.hiv:
            add_to_subsets(self.hiv_agents, agent)
            if agent.aids:
                add_to_subsets(self.hiv_aids_agents, agent)

        # Add to correct treatment set
        if agent.intervention_ever:
            add_to_subsets(self.intervention_agents, agent)
            if agent.haart:
                add_to_subsets(self.intervention_haart_agents, agent)

        if agent.prep:
            add_to_subsets(self.intervention_prep_agents, agent)

        if agent.hiv_dx:
            add_to_subsets(self.intervention_dx_agents, agent)

        if agent.incar:
            add_to_subsets(self.incarcerated_agents, agent)

        if agent.high_risk:
            add_to_subsets(self.high_risk_agents, agent)

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

    def update_agent_partners(self, graph, agent: Agent) -> bool:
        partner = get_partner(agent, self.all_agents, self.params)
        no_match = False

        def bondtype(bond_dict):
            pvalue = 0.0
            bond_probability = self.pop_random.random()
            bonded_type = "sexOnly"
            for reltype, p in bond_dict.items():
                pvalue += p
                if bond_probability < pvalue:
                    bonded_type = reltype
                    break
            return bonded_type

        if partner:
            duration = get_partnership_duration(agent, self.params)

            if agent.drug_use == "Inj" and partner.drug_use == "Inj":
                bond_type = bondtype(self.params.partnership.bond.type.PWID)
            else:
                bond_type = bondtype(self.params.partnership.bond.type[agent.so])

            relationship = Relationship(agent, partner, duration, bond_type=bond_type)

            self.relationships.append(relationship)
            graph.add_edge(
                relationship.agent1, relationship.agent2, relationship=bond_type
            )
        else:
            graph.add_node(agent)
            no_match = True

        return no_match

    def update_partner_assignments(self, graph):
        # Now create partnerships until available partnerships are out
        eligible_agents = self.all_agents
        for agent in eligible_agents:
            # add agent to network
            graph.add_node(agent)
            acquire_prob = self.params.calibration.sex.partner * (
                agent.mean_num_partners / (12.0)
            )
            if self.pop_random.random() < acquire_prob:
                self.update_agent_partners(graph, agent)
