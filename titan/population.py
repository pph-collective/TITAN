#!/usr/bin/env python
# encoding: utf-8

import random

from typing import List, Dict, Any
from scipy.stats import poisson  # type: ignore
import numpy as np  # type: ignore
from dotmap import DotMap  # type: ignore
import networkx as nx  # type: ignore

from .agent import AgentSet, Agent, Relationship
from .partnering import get_partner, get_partnership_duration
from . import utils


class Population:
    """
    :Purpose:
        This class constructs and represents the model population

    :Input:

        params : DotMap
            Model parameters

    """

    def __init__(self, params: DotMap):
        """
        :Purpose:
            Initialize Population object.
        """
        self.pop_seed = utils.get_check_rand_int(params.model.seed.ppl)

        # Init RNG for population creation to pop_seed
        self.pop_random = random.Random(self.pop_seed)
        self.np_random = np.random.RandomState(self.pop_seed)

        # this sets the global random seed for the population generation phase, during model init it gets reset at the very end
        random.seed(self.pop_seed)

        self.enable_graph = params.model.network.enable

        if self.enable_graph:
            self.graph = nx.Graph()
        else:
            self.graph = None

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
                self.add_agent(agent)

        if params.features.incar:
            self.initialize_incarceration()

        # initialize relationships
        for i in range(10):
            self.update_partner_assignments()

        if self.enable_graph:
            self.initialize_graph()

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
            sex_type = self.np_random.choice(
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
            if self.pop_random.random() < self.params.prep.pca.prep_awareness.init:
                agent.prep_awareness = True
            attprob = self.pop_random.random()
            pvalue = 0.0
            for bin, fields in self.params.pca.attitude.items():
                pvalue += fields.prob
                if attprob < pvalue:
                    agent.prep_opinion = bin
                    break

        return agent

    def add_agent(self, agent: Agent):
        """
        :Purpose:
            Create a new agent in the population.

        :Input:
            agent : int

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

        if self.enable_graph:
            self.graph.add_node(agent)

    def add_relationship(self, rel: Relationship):
        """
        :Purpose:
            Create a new relationship in the population.

        :Input:
            agent : int
        """
        self.relationships.append(rel)

        if self.enable_graph:
            self.graph.add_edge(rel.agent1, rel.agent2)

    def remove_agent(self, agent: Agent):
        """
        :Purpose:
            Remove an agent from the population.

        :Input:
            agent : int
        """
        self.all_agents.remove_agent(agent)

        if self.enable_graph:
            self.graph.remove_node(agent)

    def remove_relationship(self, rel: Relationship):
        """
        :Purpose:
            Remove a relationship from the population.

        :Input:
            agent : int
        """
        self.relationships.remove(rel)

        if self.enable_graph:
            self.graph.remove_edge(rel.agent1, rel.agent2)

    def get_age(self, race: str):
        rand = self.pop_random.random()

        bins = self.params.demographics[race].age

        for i in range(1, 6):
            if rand < bins[i].prob:
                min_age = bins[i].min
                max_age = bins[i].max
                break

        age = self.pop_random.randrange(min_age, max_age)
        return age, i

    # REVIEWED should these be in the network class? - max to incorporate with network/pop/model disentangling?

    def update_agent_partners(self, agent: Agent) -> bool:
        """
        :Purpose:
            Finds and bonds new partner. Creates relationship object for partnership, calcs
            partnership duration, and adds to networkX graph if self.enable_graph is set True.

        :Input:
            agent : Agent
            Agent that is seeking a new partner

        :Returns:
            noMatch : bool
            Bool if no match was found for agent (used for retries)
        """
        partner = get_partner(agent, self.all_agents, self.params, self.pop_random)
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
            duration = get_partnership_duration(agent, self.params, self.pop_random)

            if agent.drug_use == "Inj" and partner.drug_use == "Inj":
                bond_type = bondtype(self.params.partnership.bond.type.PWID)
            else:
                bond_type = bondtype(self.params.partnership.bond.type[agent.so])

            relationship = Relationship(agent, partner, duration, bond_type=bond_type)
            self.add_relationship(relationship)
        else:
            no_match = True

        return no_match

    def update_partner_assignments(self):
        """
        :Purpose:
            Determines which agents will seek new partners from All_agentSet.
            Calls update_agent_partners for any agents that desire partners.

        :Input:
            None
        """
        # Now create partnerships until available partnerships are out
        for agent in self.all_agents:
            acquire_prob = self.params.calibration.sex.partner * (
                agent.mean_num_partners / (12.0)
            )
            if self.pop_random.random() < acquire_prob:
                self.update_agent_partners(agent)

    def initialize_graph(self):
        """
        :Purpose:
            Initialize network with graph-based algorithm for relationship adding/pruning

        :Input:
            None
        """

        if self.params.model.network.type == "max_k_comp_size":

            def trim_component(component, max_size):
                for ag in component.nodes:
                    if self.pop_random.random() < 0.1:
                        for rel in ag.relationships:
                            if len(ag.relationships) == 1:
                                break  # Make sure that agents stay part of the network by keeping one bond
                            rel.progress(force=True)
                            self.remove_relationship(rel)
                            component.remove_edge(rel.agent1, rel.agent2)

                # recurse on new sub-components
                sub_comps = list(
                    component.subgraph(c).copy()
                    for c in nx.connected_components(component)
                )
                for sub_comp in sub_comps:
                    if sub_comp.number_of_nodes > max_size:
                        trim_component(component, max_size)

            components = sorted(self.connected_components(), key=len, reverse=True)
            for comp in components:
                if (
                    comp.number_of_nodes()
                    > self.params.model.network.component_size.max
                ):
                    print("TOO BIG", comp, comp.number_of_nodes())
                    trim_component(comp, self.params.model.network.component_size.max)

        print("Total agents in graph: ", self.graph.number_of_nodes())

    def connected_components(self):
        """
        :Purpose:
            Return connected components in graph (if enabled)

        :Input:
            agent : int
        """
        if self.enable_graph:
            return list(
                self.graph.subgraph(c).copy()
                for c in nx.connected_components(self.graph)
            )
        else:
            raise ValueError("Can't get connected_components, population doesn't have graph enabled.")
