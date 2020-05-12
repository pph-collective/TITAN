#!/usr/bin/env python
# encoding: utf-8

import random
from collections import deque
from copy import copy

from typing import List, Dict, Any, Set
import numpy as np  # type: ignore
import networkx as nx  # type: ignore
from math import ceil

from .parse_params import ObjMap
from .agent import AgentSet, Agent, Relationship
from .partnering import select_partner, get_partnership_duration
from . import utils


class Population:
    """
    :Purpose:
        This class constructs and represents the model population

    :Input:

        params : ObjMap
            Model parameters

    """

    def __init__(self, params: ObjMap):
        """
        :Purpose:
            Initialize Population object.
        """
        self.pop_seed = utils.get_check_rand_int(params.model.seed.ppl)

        # Init RNG for population creation to pop_seed
        self.pop_random = random.Random(self.pop_seed)
        self.np_random = np.random.RandomState(self.pop_seed)

        # this sets the global random seed for the population generation phase, during
        # model init it gets reset at the very end
        random.seed(self.pop_seed)

        self.enable_graph = params.model.network.enable

        if self.enable_graph:
            self.graph = nx.Graph()
        else:
            self.graph = None

        self.params = params
        # pre-fetch param sub-sets for performance
        self.demographics = params.demographics
        self.features = params.features
        self.prep = params.prep

        # build weights of population sex types by race
        self.pop_weights: Dict[str, Dict[str, List[Any]]] = {}
        self.role_weights: Dict[str, Dict] = {}
        for race in params.classes.races:
            self.role_weights[race] = {}
            self.pop_weights[race] = {}
            self.pop_weights[race]["values"] = []
            self.pop_weights[race]["weights"] = []
            for st in params.classes.sex_types:
                self.pop_weights[race]["values"].append(st)
                self.pop_weights[race]["weights"].append(
                    self.demographics[race][st].ppl
                )

                self.role_weights[race][st] = {}
                self.role_weights[race][st]["values"] = []
                self.role_weights[race][st]["weights"] = []
                for role, prob in self.demographics[race][st].sex_role.init.items():
                    self.role_weights[race][st]["values"].append(role)
                    self.role_weights[race][st]["weights"].append(prob)

        print("\tBuilding class sets")

        # All agent set list
        self.all_agents = AgentSet("AllAgents")

        # HIV status agent sets
        self.hiv_agents = AgentSet("HIV", parent=self.all_agents)

        # High risk agent sets
        self.high_risk_agents = AgentSet("HRisk", parent=self.all_agents)

        # pwid agents (performance for partnering)
        self.pwid_agents = AgentSet("PWID", parent=self.all_agents)

        # agents who can take on a partner
        self.partnerable_agents: Dict[str, Set[Agent]] = {}
        for bond_type in self.params.classes.bond_types.keys():
            self.partnerable_agents[bond_type] = set()

        # who can sleep with whom
        self.sex_partners: Dict[str, Set[Agent]] = {}
        for sex_type in self.params.classes.sex_types.keys():
            self.sex_partners[sex_type] = set()

        self.relationships: Set[Relationship] = set()

        # keep track of prep agent counts by race
        self.prep_counts = {race: 0 for race in params.classes.races}

        # find average partnership durations
        self.mean_rel_duration: Dict[str, int] = {}
        for bond in self.params.partnership.duration:
            if self.params.partnership.duration[bond].type == "bins":
                weights = []
                vals = []
                dur_bins = params.partnership.duration[bond].bins
                for bins in dur_bins:
                    if bins > 1:
                        weights.append(dur_bins[bins].prob - dur_bins[bins - 1].prob)
                    else:
                        weights.append(dur_bins[bins].prob)
                    vals.append(np.average([dur_bins[bins].min, dur_bins[bins].max]))
                self.mean_rel_duration[bond] = np.average(vals, weights=weights)
            else:
                self.mean_rel_duration[bond] = self.params.partnership.duration[
                    bond
                ].distribution.mean

        print("\tCreating agents")

        for race in params.classes.races:
            for i in range(round(params.model.num_pop * params.demographics[race].ppl)):
                agent = self.create_agent(race)
                self.add_agent(agent)

        if params.features.incar:
            print("\tInitializing Incarceration")
            self.initialize_incarceration()

        # initialize relationships
        print("\tCreating Relationships")
        self.update_partner_assignments()

        if self.enable_graph:
            self.initialize_graph()

    def initialize_incarceration(self):

        for a in self.all_agents.members:
            incar_params = self.demographics[a.race][a.so].incar
            jail_duration = incar_params.duration.init

            prob_incar = incar_params.init
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
        if self.pop_random.random() < self.demographics[race]["PWID"].ppl:
            drug_type = "Inj"
        else:
            drug_type = "None"

        age, age_bin = self.get_age(race)

        agent = Agent(sex_type, age, race, drug_type)
        agent.age_bin = age_bin
        agent.sex_role = self.np_random.choice(
            self.role_weights[race][sex_type]["values"],
            p=self.role_weights[race][sex_type]["weights"],
        )
        if self.features.msmw and sex_type == "HM":
            if self.pop_random.random() < self.params.msmw.prob:
                agent.msmw = True

        if drug_type == "Inj":
            agent_params = self.demographics[race]["PWID"]
        else:
            agent_params = self.demographics[race][sex_type]

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

                    haart_adh = self.demographics[race][sex_type].haart.adherence
                    if self.pop_random.random() < haart_adh:
                        adherence = 5
                    else:
                        adherence = self.pop_random.randint(1, 4)

                    # add to agent haart set
                    agent.haart_adherence = adherence
                    agent.haart_time = 0

            # if HIV, how long has the agent had it? Random sample
            agent.hiv_time = self.pop_random.randint(1, self.params.hiv.max_init_time)

        else:
            if self.features.prep:
                if self.prep.start == 0:
                    prob_prep = self.prep.target
                else:
                    prob_prep = 0.0

                if self.pop_random.random() < prob_prep:
                    agent.prep = True
                    agent.intervention_ever = True
                    if (
                        self.pop_random.random() > self.prep.lai.prob
                        and "Inj" in self.prep.type
                    ):
                        agent.prep_type = "Inj"
                    else:
                        agent.prep_type = "Oral"

        # Check if agent is HR as baseline.
        if (
            self.features.high_risk
            and self.pop_random.random()
            < self.demographics[race][sex_type].high_risk.init
        ):
            agent.high_risk = True
            agent.high_risk_ever = True
            agent.high_risk_time = self.pop_random.randint(
                1, self.params.high_risk.sex_based[agent.so].duration
            )

        # get agent's mean partner numbers for bond type
        def partner_distribution(dist_info):
            return ceil(
                utils.safe_dist(dist_info, self.np_random)
                * utils.safe_divide(
                    self.params.calibration.sex.partner, self.mean_rel_duration[bond]
                )
            )

        for bond in self.params.classes.bond_types.keys():
            agent.partners[bond] = set()
            if agent.drug_use == "Inj":
                dist_info = self.params.demographics[agent.race].PWID.num_partners[bond]
            else:
                dist_info = self.params.demographics[agent.race][agent.so].num_partners[
                    bond
                ]
            agent.mean_num_partners[bond] = partner_distribution(dist_info)
            # so not zero if added mid-year
            agent.target_partners[bond] = agent.mean_num_partners[bond]
            if bond == "Inj":
                assert agent.drug_use == "Inj" or agent.mean_num_partners[bond] == 0

            if agent.target_partners[bond] > 0:
                self.partnerable_agents[bond].add(agent)

        if self.features.pca:
            if self.pop_random.random() < self.prep.pca.awareness.init:
                agent.prep_awareness = True
            attprob = self.pop_random.random()
            pvalue = 0.0
            for bin, fields in self.prep.pca.attitude.items():
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

        # Add to all agent set
        self.all_agents.add_agent(agent)

        if agent.hiv:
            self.hiv_agents.add_agent(agent)

        if agent.high_risk:
            self.high_risk_agents.add_agent(agent)

        if agent.drug_use == "Inj":
            self.pwid_agents.add_agent(agent)

        # who can sleep with this agent
        for sex_type in self.params.classes.sex_types[agent.so].sleeps_with:
            self.sex_partners[sex_type].add(agent)

        if agent.prep:
            self.prep_counts[agent.race] += 1

        if self.enable_graph:
            self.graph.add_node(agent)

    def add_relationship(self, rel: Relationship):
        """
        :Purpose:
            Create a new relationship in the population.

        :Input:
            agent : int
        """
        self.relationships.add(rel)

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

        for partner_type in self.sex_partners:
            if agent in self.sex_partners[partner_type]:
                self.sex_partners[partner_type].remove(agent)

        if agent.prep:
            self.prep_counts[agent.race] -= 1

        if self.enable_graph:
            self.graph.remove_node(agent)

        for bond in self.partnerable_agents.values():
            if agent in bond:
                bond.remove(agent)

    def remove_relationship(self, rel: Relationship):
        """
        :Purpose:
            Remove a relationship from the population.

        :Input:
            agent : int
        """
        self.relationships.remove(rel)

        # without this relationship, are agents partnerable again?
        self.update_partnerability(rel.agent1)
        self.update_partnerability(rel.agent2)

        if self.enable_graph:
            self.graph.remove_edge(rel.agent1, rel.agent2)

    def get_age(self, race: str):
        """
        :Purpose:
            Get an age of an agent, given their race

        :Input:
            race : str

        :Returns:
            age : int
            bin : int
        """
        rand = self.pop_random.random()

        bins = self.demographics[race].age

        for i in range(1, 6):
            if rand < bins[i].prob:
                min_age = bins[i].min
                max_age = bins[i].max
                break

        age = self.pop_random.randrange(min_age, max_age)
        return age, i

    def update_agent_partners(self, agent: Agent, bond_type: str) -> bool:
        """
        :Purpose:
            Finds and bonds new partner. Creates relationship object for partnership,
            calcs partnership duration, and adds to networkX graph if self.enable_graph
            is set True.

        :Input:
            agent : Agent
            Agent that is seeking a new partner

        :Returns:
            noMatch : bool
            Bool if no match was found for agent (used for retries)
        """
        partner = select_partner(
            agent,
            self.partnerable_agents[bond_type],
            self.sex_partners,
            self.pwid_agents,
            self.params,
            self.pop_random,
            bond_type,
        )
        no_match = True

        if partner:
            duration = get_partnership_duration(self.params, self.np_random, bond_type)
            relationship = Relationship(agent, partner, duration, bond_type=bond_type)
            self.add_relationship(relationship)
            # can partner still partner?
            if len(partner.partners[bond_type]) >= (
                partner.target_partners[bond_type]
                * self.params.calibration.partnership.buffer
            ):
                self.partnerable_agents[bond_type].remove(partner)
            no_match = False
        return no_match

    def update_partner_assignments(self, t=0):
        """
        :Purpose:
            Determines which agents will seek new partners from All_agentSet.
            Calls update_agent_partners for any agents that desire partners.

        :Input:
            None
        """
        print(f"target partnerships: ")

        # update agent targets annually
        if t % self.params.model.time.steps_per_year == 0:
            self.update_partner_targets()

        # Now create partnerships until available partnerships are out
        for bond in self.params.classes.bond_types:
            eligible_agents = deque(
                [
                    a
                    for a in self.all_agents
                    if len(a.partners[bond]) < a.target_partners[bond]
                ]
            )
            attempts = {a: 0 for a in eligible_agents}

            while eligible_agents:
                agent = eligible_agents.popleft()

                # no match
                if self.update_agent_partners(agent, bond):
                    attempts[agent] += 1

                # add agent back to eligible pool
                if (
                    len(agent.partners[bond]) < agent.target_partners[bond]
                    and attempts[agent]
                    < self.params.calibration.partnership.break_point
                ):
                    eligible_agents.append(agent)

        print(
            f"actual partnerships (post): "
            f"{sum([len(a.partners) for a in self.all_agents])}"
        )

    def update_partner_targets(self):
        for a in self.all_agents:
            for bond in self.params.classes.bond_types:
                a.target_partners[bond] = utils.poisson(
                    a.mean_num_partners[bond], self.np_random
                )
            self.update_partnerability(a)

    def update_partnerability(self, a):
        # update partnerability
        for bond in self.params.classes.bond_types.keys():
            if a in self.partnerable_agents[bond]:
                if len(a.partners[bond]) > (
                    a.target_partners[bond] * self.params.calibration.partnership.buffer
                ):
                    self.partnerable_agents[bond].remove(a)
            elif len(a.partners) < (
                a.target_partners[bond] * self.params.calibration.partnership.buffer
            ):
                self.partnerable_agents[bond].add(a)

    def initialize_graph(self):
        """
        :Purpose:
            Initialize network with graph-based algorithm for relationship
            adding/pruning

        :Input:
            None
        """

        if self.params.model.network.type == "max_k_comp_size":

            def trim_component(component, max_size):
                for ag in component.nodes:
                    if (
                        self.pop_random.random()
                        < self.params.calibration.network.trim.prob
                    ):
                        for rel in copy(ag.relationships):
                            if len(ag.relationships) == 1:
                                break  # Make sure that agents stay part of the
                                # network by keeping one bond
                            rel.progress(force=True)
                            self.remove_relationship(rel)
                            component.remove_edge(rel.agent1, rel.agent2)

                # recurse on new sub-components
                sub_comps = list(
                    component.subgraph(c).copy()
                    for c in nx.connected_components(component)
                )
                for sub_comp in sub_comps:
                    if sub_comp.number_of_nodes() > max_size:
                        trim_component(component, max_size)

            components = sorted(self.connected_components(), key=len, reverse=True)
            for comp in components:
                if (
                    comp.number_of_nodes()
                    > self.params.model.network.component_size.max
                ):
                    print("TOO BIG", comp, comp.number_of_nodes())
                    trim_component(comp, self.params.model.network.component_size.max)

        print("\tTotal agents in graph: ", self.graph.number_of_nodes())

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
            raise ValueError(
                "Can't get connected_components, population doesn't have graph enabled."
            )
