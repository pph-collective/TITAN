#!/usr/bin/env python
# encoding: utf-8

import random
from collections import deque
from copy import copy
from math import ceil
from typing import List, Dict, Set, Optional, Tuple
import logging

import numpy as np  # type: ignore
import networkx as nx  # type: ignore
import nanoid  # type: ignore

from . import parse_params
from . import agent as ag
from . import location
from . import partnering
from . import utils
from . import features
from . import exposures


class Population:
    def __init__(self, params: "parse_params.ObjMap", id: Optional[str] = None):
        """
        Initialize Population object.

        args:
            params : Model parameters
            id: 8 character identifier for a model
        """
        if id is None:
            self.id = nanoid.generate(size=8)
        else:
            self.id = id

        utils.set_up_logging(params)

        logging.info(f"Population ID: {self.id}")

        self.pop_seed = utils.get_check_rand_int(params.model.seed.ppl)

        # Init RNG for population creation to pop_seed
        self.pop_random = random.Random(self.pop_seed)
        self.np_random = np.random.default_rng(self.pop_seed)

        self.enable_graph = params.model.network.enable
        self.components: List = []

        if self.enable_graph:
            self.graph = nx.Graph()
        else:
            self.graph = None

        self.params = params

        # set up the in-scope exposures
        self.exposures = [
            exposure
            for exposure in exposures.BaseExposure.__subclasses__()
            if self.params.exposures[exposure.name]
        ]
        # initialize the class level items
        for exposure in self.exposures:
            exposure.init_class(params)

        # set up the in-scope features
        self.features = [
            feature
            for feature in features.BaseFeature.__subclasses__()
            if self.params.features[feature.name]
        ]
        # initialize the class level items
        for feature in self.features:
            feature.init_class(params)

        # set up the population's locations and edges
        self.geography = location.Geography(params)

        # All agent set list
        self.all_agents = ag.AgentSet("AllAgents")

        # pwid agents (performance for partnering)
        self.pwid_agents = ag.AgentSet("PWID", parent=self.all_agents)

        # agents who can take on a partner
        self.partnerable_agents: Dict[str, Set["ag.Agent"]] = {}
        for bond_type in self.params.classes.bond_types.keys():
            self.partnerable_agents[bond_type] = set()

        # who can sleep with whom
        self.sex_partners: Dict[str, Set["ag.Agent"]] = {}
        for sex_type in self.params.classes.sex_types.keys():
            self.sex_partners[sex_type] = set()

        self.relationships: Set["ag.Relationship"] = set()

        # find average partnership durations
        self.mean_rel_duration: Dict[str, Dict] = partnering.get_mean_rel_duration(
            self.params
        )

        logging.info("  Creating agents")
        # for each location in the population, create agents per that location's demographics
        init_time = -1 * self.params.model.time.burn_steps
        for loc in self.geography.locations.values():
            for race in params.classes.races:
                for i in range(
                    round(
                        params.model.num_pop
                        * loc.ppl
                        * loc.params.demographics[race].ppl
                    )
                ):
                    if self.all_agents.num_members() >= self.params.model.num_pop:
                        logging.warning(
                            "WARNING: not adding agent to population - too many agents"
                        )
                        break
                    agent = self.create_agent(loc, race, init_time)
                    self.add_agent(agent)

        # initialize relationships
        logging.info("  Creating Relationships")
        self.update_partner_assignments(0)

    def create_agent(
        self,
        loc: "location.Location",
        race: str,
        time: int,
        sex_type: Optional[str] = None,
        drug_type: Optional[str] = None,
        age: Optional[int] = None,
    ) -> "ag.Agent":
        """
        Create a new agent with randomly assigned attributes according to population
        demographics [params.demographics]

        args:
            loc: location the agent will live in
            race: race of the new agent
            time: current time step of the model
            sex_type: sex_type of the new agent
            drug_type: drug_type of the new agent
            age: age of the new agent

        returns:
             a new agent
        """

        if not sex_type:
            sex_type = utils.safe_random_choice(
                loc.pop_weights[race]["values"],
                self.pop_random,
                weights=loc.pop_weights[race]["weights"],
            )
            # no choice available
            if sex_type is None:
                raise ValueError("Agent must have sex type")

        # Determine drugtype
        if not drug_type:
            drug_type = utils.safe_random_choice(
                loc.drug_weights[race][sex_type]["values"],
                self.pop_random,
                weights=loc.drug_weights[race][sex_type]["weights"],
            )
            # no choice available
            if drug_type is None:
                raise ValueError("Agent must have drug type")

        if not age:
            age = self.get_age(loc, race)

        agent = ag.Agent(sex_type, age, race, drug_type, loc)

        sex_role = utils.safe_random_choice(
            loc.role_weights[race][sex_type]["values"],
            self.pop_random,
            weights=loc.role_weights[race][sex_type]["weights"],
        )
        if sex_role is None:
            raise ValueError("Agent must have sex role")
        else:
            agent.sex_role = sex_role

        agent_params = (
            agent.location.params.demographics[race]
            .sex_type[sex_type]
            .drug_type[drug_type]
        )

        for exposure in self.exposures:
            agent_feature = getattr(agent, exposure.name)
            agent_feature.init_agent(self, time)

        for bond, bond_def in loc.params.classes.bond_types.items():
            agent.partners[bond] = set()
            dist_info = agent_params.num_partners[bond]
            agent.mean_num_partners[bond] = ceil(
                utils.safe_dist(dist_info, self.np_random)
                * utils.safe_divide(
                    agent.location.params.calibration.sex.partner,
                    self.mean_rel_duration[bond][race],
                )
            )
            # so not zero if added mid-year
            agent.target_partners[bond] = agent.mean_num_partners[bond]
            if "injection" in bond_def.acts_allowed:
                assert agent.drug_type == "Inj" or agent.mean_num_partners[bond] == 0

            if agent.target_partners[bond] > 0:
                self.partnerable_agents[bond].add(agent)

        for feature in self.features:
            agent_feature = getattr(agent, feature.name)
            agent_feature.init_agent(self, time)

        return agent

    def add_agent(self, agent: "ag.Agent"):
        """
        Adds an agent to the population

        args:
            agent : The agent to be added
        """
        # Add to all agent set
        self.all_agents.add_agent(agent)

        if agent.drug_type == "Inj":
            self.pwid_agents.add_agent(agent)

        # who can sleep with this agent
        for sex_type in self.params.classes.sex_types[agent.sex_type].sleeps_with:
            self.sex_partners[sex_type].add(agent)

        if self.enable_graph:
            self.graph.add_node(agent)

    def add_relationship(self, rel: "ag.Relationship"):
        """
        Add a new relationship to the population.

        args:
            rel : The Relationship to be added
        """
        self.relationships.add(rel)

        if self.enable_graph:
            self.graph.add_edge(rel.agent1, rel.agent2, type=rel.bond_type)

    def remove_agent(self, agent: "ag.Agent"):
        """
        Remove an agent from the population.

        args:
            agent : Agent to remove
        """
        for rel in copy(agent.relationships):
            rel.progress(force=True)
            self.remove_relationship(rel)

        self.all_agents.remove_agent(agent)

        for partner_type in self.sex_partners:
            if agent in self.sex_partners[partner_type]:
                self.sex_partners[partner_type].remove(agent)

        for exposure in self.exposures:
            agent_attr = getattr(agent, exposure.name)
            if agent_attr.active:
                exposure.remove_agent(agent)

        for feature in self.features:
            agent_attr = getattr(agent, feature.name)
            if agent_attr.active:
                feature.remove_agent(agent)

        if self.enable_graph:
            self.graph.remove_node(agent)

        for bond in self.partnerable_agents.values():
            if agent in bond:
                bond.remove(agent)

        # mark agent component as -1 (no component)
        agent.component = "-1"

    def remove_relationship(self, rel: "ag.Relationship"):
        """
        Remove a relationship from the population.

        args:
            rel : Relationship to remove
        """
        self.relationships.remove(rel)

        # without this relationship, are agents partnerable again?
        self.update_partnerability(rel.agent1)
        self.update_partnerability(rel.agent2)

        if self.enable_graph:
            self.graph.remove_edge(rel.agent1, rel.agent2)

    def get_age(self, loc: "location.Location", race: str) -> int:
        """
        Given the population characteristics, get a random age to assign to an agent given the race of that agent

        args:
            race : race of the agent whose age is being generated

        returns:
            age and the bin the age came from
        """
        bins = loc.params.demographics[race].age
        i = utils.get_independent_bin(self.pop_random, bins)
        age = self.pop_random.randrange(bins[i].min, bins[i].max)
        return age

    def update_agent_partners(
        self, agent: "ag.Agent", bond_type: str, components: List
    ) -> bool:
        """
        Finds and bonds new partner. Creates relationship object for partnership,
            calcs partnership duration, adds it to the population, and adds to networkX graph if self.enable_graph
            is set True.

        args:
            agent: Agent that is seeking a new partner
            bond_type: What type of bond the agent is seeking to make

        returns:
            True if no match was found for agent (used for retries)
        """
        partnerable_agents = self.partnerable_agents[bond_type]
        if (
            self.pop_random.random()
            < self.params.partnership.network.same_component.prob
            and agent.has_partners()
        ):
            # find agent's component
            agent_component: Set["ag.Agent"] = set()
            for comp in components:
                if agent in comp:
                    agent_component = comp
                    break

            partnerable_agents = partnerable_agents & agent_component

        partner = partnering.select_partner(
            agent,
            partnerable_agents,
            self.sex_partners,
            self.pwid_agents,
            self.params,
            self.pop_random,
            bond_type,
        )
        no_match = True

        if partner:
            race = utils.safe_random_choice([agent.race, partner.race], self.pop_random)
            duration = partnering.get_partnership_duration(
                agent.location.params, self.np_random, bond_type, race
            )
            relationship = ag.Relationship(
                agent, partner, duration, bond_type=bond_type
            )
            self.add_relationship(relationship)
            # can partner still partner?
            if len(partner.partners[bond_type]) > (
                partner.target_partners[bond_type]
                * self.params.calibration.partnership.buffer
            ):
                self.partnerable_agents[bond_type].remove(partner)
            no_match = False
        return no_match

    def update_partner_assignments(self, t: int):
        """
        Determines which agents will seek new partners from All_agentSet.
            Calls update_agent_partners for any agents that desire partners.

        args:
            t: current time step of the model
        """
        # update agent targets annually
        if t % self.params.model.time.steps_per_year == 0:
            self.update_partner_targets()

        if self.enable_graph:
            network_components = [set(g.nodes()) for g in self.components]
        else:
            network_components = []

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
                if len(agent.partners[bond]) < agent.target_partners[bond]:

                    # no match
                    if self.update_agent_partners(agent, bond, network_components):
                        attempts[agent] += 1

                    # add agent back to eligible pool
                    if (
                        len(agent.partners[bond]) < agent.target_partners[bond]
                        and attempts[agent]
                        < self.params.calibration.partnership.break_point
                    ):
                        eligible_agents.append(agent)

        if self.enable_graph:
            self.trim_graph()

        self.update_agent_components()

    def update_partner_targets(self):
        """
        Update the target number of partners for each agent and bond type
        """
        for a in self.all_agents:
            for bond in self.params.classes.bond_types:
                a.target_partners[bond] = utils.poisson(
                    a.mean_num_partners[bond], self.np_random
                )
            self.update_partnerability(a)

    def update_partnerability(self, a):
        """
        Update whether each agent in the population is currently able to form new relationships for each bond type
        """
        for bond in self.params.classes.bond_types.keys():
            if a in self.partnerable_agents[bond]:
                if len(a.partners[bond]) > (
                    a.target_partners[bond] * self.params.calibration.partnership.buffer
                ):
                    self.partnerable_agents[bond].remove(a)
            elif len(a.partners[bond]) < (
                a.target_partners[bond] * self.params.calibration.partnership.buffer
            ):
                self.partnerable_agents[bond].add(a)

    def update_agent_components(self):
        """
        Update the component IDs associated with each agent based on the current state of the graph
        """
        if self.enable_graph:
            self.components = utils.connected_components(self.graph)
            for id, component in enumerate(self.components):
                for agent in component.nodes:
                    agent.component = str(id)

            self.params.classes.components = list(
                map(str, range(-1, len(self.components)))
            )

    def trim_graph(self):
        """
        Initialize network with graph-based algorithm for relationship
            adding/pruning
        """
        if self.params.model.network.type == "comp_size":

            def trim_component(component, max_size):
                for agent in component.nodes:
                    if (
                        self.pop_random.random()
                        < self.params.calibration.network.trim.prob
                    ):
                        for rel in copy(agent.relationships):
                            if len(agent.relationships) == 1:
                                break  # Make sure that agents stay part of the
                                # network by keeping one bond
                            rel.progress(force=True)
                            self.remove_relationship(rel)

                # recurse on new sub-components
                sub_comps = utils.connected_components(component)
                for sub_comp in sub_comps:
                    if sub_comp.number_of_nodes() > max_size:
                        trim_component(component, max_size)
                    else:
                        break

            components = self.connected_components()
            for comp in components:
                if (
                    comp.number_of_nodes()
                    > self.params.model.network.component_size.max
                ):
                    logging.info("TOO BIG", comp, comp.number_of_nodes())
                    trim_component(comp, self.params.model.network.component_size.max)

        logging.info(f"  Total agents in graph: {self.graph.number_of_nodes()}")

    def connected_components(self) -> List:
        """
        Get connected components in graph (if enabled)

        returns:
            list of connected components
        """
        if self.enable_graph:
            return self.components
        else:
            raise ValueError(
                "Can't get connected_components, population doesn't have graph enabled."
            )

    def migrate(self):
        """
        Have agents migrate between locations with probabilities defined in `location.migration.matrix_file`.
        """
        m_attr = self.params.location.migration.attribute
        for a in self.all_agents:
            m_param = a.location.migration_weights
            if self.pop_random.random() < m_param["prob"]:
                new_loc = utils.safe_random_choice(
                    m_param["values"],
                    self.pop_random,
                    weights=m_param["weights"],
                )
                if m_attr == "name":
                    a.location = self.geography.locations[new_loc]
                elif m_attr == "category":
                    a.location = utils.safe_random_choice(
                        self.geography.categories[new_loc], self.pop_random
                    )
