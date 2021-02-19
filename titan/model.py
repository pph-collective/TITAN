import random
from typing import Dict, List, Optional
from copy import copy
import os

import numpy as np  # type: ignore
import nanoid  # type: ignore

from . import agent as ag
from . import output as ao
from . import probabilities as prob
from .parse_params import ObjMap
from . import exposures, features, interactions, population, utils


class TITAN:
    def __repr__(self):
        res = "\n"
        res += f"Seed: {self.run_seed}\n"
        res += f"Npop: {self.params.model.num_pop}\n"
        res += f"Time: {self.params.model.time.num_steps}\n"

        return res

    def __init__(
        self,
        params: ObjMap,
        pop: Optional["population.Population"] = None,
    ):
        """
        This is the core class used to simulate the spread of exposures through a relationship based network.

        args:
            params: the parameter object for this model
            pop: an initialized population to run the model on
        """

        self.params = params
        # pre-fetch commonly used param sub-sets for performance
        self.calibration = params.calibration

        print("=== Begin Initialization Protocol ===\n")

        if pop is None:
            print("\tGenerating new population")
            self.pop = population.Population(params)
        else:
            print("\tUsing provided population")
            self.pop = pop

        self.time = -1 * self.params.model.time.burn_steps  # burn is negative time
        self.id = nanoid.generate(size=8)

        self.features = [
            feature
            for feature in features.BaseFeature.__subclasses__()
            if self.params.features[feature.name]
        ]

        # set up the in-scope exposures
        self.exposures = [
            exposure
            for exposure in exposures.BaseExposure.__subclasses__()
            if self.params.exposures[exposure.name]
        ]

        self.interactions = {
            interaction.name: interaction
            for interaction in interactions.BaseInteraction.__subclasses__()
        }

        # Set seed format. 0: pure random, else: fixed value
        self.run_seed = utils.get_check_rand_int(params.model.seed.run)
        print(f"\tRun seed was set to: {self.run_seed}")
        self.run_random = random.Random(self.run_seed)
        self.np_random = np.random.default_rng(self.run_seed)
        random.seed(self.run_seed)
        print(("\tFIRST RANDOM CALL {}".format(random.randint(0, 100))))

        print("\tResetting death count")
        self.deaths: List["ag.Agent"] = []  # Number of death

        print("\n === Initialization Protocol Finished ===")

    def print_stats(self, stat: Dict[str, Dict[str, int]], outdir: str):
        """
        Create/update all of the reports defined in the params
        """
        for report in self.params.outputs.reports:
            printer = getattr(ao, report)
            printer(
                self.id,
                self.time,
                self.run_seed,
                self.pop.pop_seed,
                stat,
                self.params,
                outdir,
            )

        # network-based reports
        if (
            self.time % self.params.outputs.print_frequency == 0
            and self.params.model.network.enable
        ):
            network_outdir = os.path.join(outdir, "network")
            if self.params.outputs.network.calc_component_stats:
                ao.print_components(
                    self.id,
                    self.time,
                    self.run_seed,
                    self.pop.pop_seed,
                    self.pop.connected_components(),
                    network_outdir,
                )

            if self.params.outputs.network.calc_network_stats:
                ao.write_network_stats(
                    self.pop.graph, network_outdir, self.id, self.time
                )

            if self.params.outputs.network.edge_list:
                ao.write_graph_edgelist(
                    self.pop.graph, network_outdir, self.id, self.time
                )

    def reset_trackers(self):
        self.deaths = []

    def run(self, outdir: str):
        """
        Runs the model for the number of time steps defined in params, at each time step does:

        1. Increments time
        2. Takes one step
        3. Resets trackers

        args:
            outdir: path to directory where results should be saved
        """
        # make sure initial state of things get printed
        stats = ao.get_stats(
            self.pop.all_agents,
            self.deaths,
            self.params,
            self.exposures,
            self.features,
            self.time,
        )
        self.print_stats(stats, outdir)

        if self.params.model.time.burn_steps > 0:
            print("\t===! Start Burn Loop !===")

        # time starts at negative burn steps, model run starts at t = 1
        while self.time < self.params.model.time.num_steps:
            if self.time == 0:
                if self.params.model.time.burn_steps > 0:
                    print("\t===! Burn Loop Complete !===")
                print("\t===! Start Main Loop !===")

            self.time += 1
            self.step(outdir)
            self.reset_trackers()

        print("\t===! Main Loop Complete !===")

    def step(self, outdir: str):
        """
        A single time step in the model:

        1. Perform timeline_scaling updates to params if needed
        2. Update all agents
        3. Write/update reports with this timestep's data

        args:
            outdir: path to directory where reports should be saved
        """
        print(f"\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t.: TIME {self.time}")
        print(
            "\tSTARTING HIV count:{}\tTotal Incarcerated:{}\tHR+:{}\t"
            "PrEP:{}".format(
                len(exposures.HIV.agents),
                sum([1 for a in self.pop.all_agents if a.incar.active]),  # type: ignore[attr-defined]
                sum([1 for a in self.pop.all_agents if a.high_risk.active]),  # type: ignore[attr-defined]
                sum([1 for a in self.pop.all_agents if a.prep.active]),  # type: ignore[attr-defined]
            )
        )

        self.timeline_scaling()

        self.update_all_agents()

        stats = ao.get_stats(
            self.pop.all_agents,
            self.deaths,
            self.params,
            self.exposures,
            self.features,
            self.time,
        )
        self.print_stats(stats, outdir)

        print(f"Number of relationships: {len(self.pop.relationships)}")
        self.pop.all_agents.print_subsets()

    def update_all_agents(self):
        """
        The core of the model.  For a time step, update all of the agents and relationships:


        1. End relationships with no remaining duration
        2. Agent death/replacement
        3. Update partner assignments (create new relationships as needed)
        4. Create an agent zero (if enabled and the time is right)
        5. Agents in relationships interact
        6. Update features at the population level
        7. Update each agent's status for:
            * age
            * all exposures
            * all features (agent level)
        """
        # If static network, ignore relationship progression
        if not self.params.features.static_network:
            for rel in copy(self.pop.relationships):
                if rel.progress():
                    self.pop.remove_relationship(rel)

        if self.params.features.die_and_replace:
            self.die_and_replace()

        if not self.params.features.static_network:
            self.pop.update_partner_assignments(t=self.time)

        # If agent zero enabled, create agent zero at the beginning of main loop.
        if (
            self.time == self.params.agent_zero.start_time
            and self.params.features.agent_zero
        ):
            self.make_agent_zero()

        for rel in self.pop.relationships:
            self.agents_interact(rel)

        for feature in self.features:
            feature.update_pop(self)

        for agent in self.pop.all_agents:
            self.update_agent(agent)

    def update_agent(self, agent):
        """
        Update an agent at the given model timestep.

        Update the agent's status for:
            * age
            * all exposures
            * all features (agent level)
        """
        # happy birthday agents!
        if self.time > 0 and (self.time % self.params.model.time.steps_per_year) == 0:
            agent.age += 1

        for exposure in self.exposures:
            agent_feature = getattr(agent, exposure.name)
            agent_feature.update_agent(self)

        for feature in self.features:
            agent_feature = getattr(agent, feature.name)
            agent_feature.update_agent(self)

    def make_agent_zero(self):
        """
        Identify an agent as agent zero and HIV convert them
        """
        bonds = [  # Find what bond_types have the allowed interaction
            bond
            for bond, act_type in self.params.classes.bond_types.items()
            if self.params.agent_zero.interaction_type in act_type.acts_allowed
        ]
        max_partners = 0
        max_agent = None
        zero_eligible = []
        for agent in self.pop.all_agents:
            num_partners = agent.get_num_partners(bond_types=bonds)
            if num_partners >= self.params.agent_zero.num_partners:
                zero_eligible.append(agent)
            if num_partners > max_partners:
                max_partners = num_partners
                max_agent = agent

        agent_zero = utils.safe_random_choice(zero_eligible, self.run_random)
        if agent_zero:  # if eligible agent, make agent 0
            zero_attr = getattr(agent_zero, self.params.agent_zero.exposure)
            zero_attr.convert(self)
        elif self.params.agent_zero.fallback and max_agent is not None:
            zero_attr = getattr(max_agent, self.params.agent_zero.exposure)
            zero_attr.convert(self)
        else:
            raise ValueError("No agent zero!")

    def timeline_scaling(self):
        """
        Scale/un-scale any params with timeline_scaling definitions per their
        definition.  Applied to all parameters (main model, and location specific).
        """
        if not self.params.features.timeline_scaling:
            return None

        # gather all of the param objectss to be scaled
        params_set = {self.params}
        for location in self.pop.geography.locations.values():
            params_set.add(location.params)

        # iterate over each param and update the values if the time is right
        for params in params_set:
            for defn in params.timeline_scaling.timeline.values():
                param = defn.parameter
                if param != "ts_default":
                    if defn.start_time == self.time:
                        print(f"timeline scaling - {param}")
                        utils.scale_param(params, param, defn.scalar)
                    elif defn.stop_time == self.time:
                        print(f"timeline un-scaling - {param}")
                        utils.scale_param(params, param, 1 / defn.scalar)

    def agents_interact(self, rel: "ag.Relationship"):
        """
        Let an agent interact with a partner.

        Based on the interaction types of the relationship, interact in the following ways:

        * Peer Change Agent
        * Injection
        * Sex

        args:
            rel : The relationship that the agents interact in
        """
        interaction_types = self.params.classes.bond_types[rel.bond_type].acts_allowed
        # If either agent is incarcerated, skip their interaction
        if rel.agent1.incar.active or rel.agent2.incar.active:  # type: ignore[attr-defined]
            return

        for interaction_type in interaction_types:
            interaction = self.interactions[interaction_type]
            interaction.interact(self, rel)

    def die_and_replace(self):
        """
        Let agents die and replace the dead agent with a new agent randomly.
        """
        # die stage
        for agent in self.pop.all_agents:
            # agent incarcerated, don't evaluate for death
            if agent.incar.active:
                continue

            # death rate per 1 person-month
            p = (
                prob.get_death_rate(
                    agent.hiv.active,
                    agent.hiv.aids,
                    agent.drug_type,
                    agent.sex_type,
                    agent.haart.adherent,
                    agent.race,
                    agent.location,
                    self.params.model.time.steps_per_year,
                )
                * self.calibration.mortality
            )

            if self.run_random.random() < p:
                self.deaths.append(agent)

                # End all existing relationships
                for rel in copy(agent.relationships):
                    rel.progress(force=True)
                    self.pop.remove_relationship(rel)

        # replace stage
        for agent in self.deaths:
            # mark agent component as -1 (no componenet)
            agent.component = "-1"

            # Remove agent from agent class and sub-sets
            self.pop.remove_agent(agent)

            new_agent = self.pop.create_agent(
                agent.location, agent.race, self.time, agent.sex_type, agent.drug_type
            )
            self.pop.add_agent(new_agent)
