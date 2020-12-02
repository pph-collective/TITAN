# Imports
import random
from typing import Dict, List, Optional
from copy import copy
import os

import numpy as np  # type: ignore
import nanoid  # type: ignore


from . import agent as ag
from . import population
from .network import NetworkGraphUtils
from . import output as ao
from . import probabilities as prob
from . import utils
from .parse_params import ObjMap
from . import features
from . import interactions


class HIVModel:
    def __repr__(self):
        res = "\n"
        res += f"Seed: {self.run_seed}\n"
        res += f"Npop: {self.params.model.num_pop}\n"
        res += f"Time: {self.params.model.time.num_steps}\n"

        return res

    def __init__(
        self, params: ObjMap, pop: Optional["population.Population"] = None,
    ):
        """
        This is the core class used to simulate
            the spread of HIV and drug use in one geography.

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

        self.network_utils: Optional[NetworkGraphUtils]
        if params.model.network.enable:
            self.network_utils = NetworkGraphUtils(self.pop.graph)
        else:
            self.network_utils = None

        self.time = -1 * self.params.model.time.burn_steps  # burn is negative time
        self.id = nanoid.generate(size=8)

        self.features = [
            feature
            for feature in features.BaseFeature.__subclasses__()
            if self.params.features[feature.name]
        ]

        self.interactions = {
            interaction.name: interaction
            for interaction in interactions.BaseInteraction.__subclasses__()
        }

        # Set seed format. 0: pure random, else: fixed value
        self.run_seed = utils.get_check_rand_int(params.model.seed.run)
        print(f"\tRun seed was set to: {self.run_seed}")
        self.run_random = random.Random(self.run_seed)
        self.np_random = np.random.RandomState(self.run_seed)
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
            assert (
                self.network_utils is not None
            ), "Graph must be enabled to print network reports"

            network_outdir = os.path.join(outdir, "network")
            if self.params.outputs.network.draw_figures:
                self.network_utils.visualize_network(
                    network_outdir, curtime=self.time, label=f"{self.id}"
                )

            if self.params.outputs.network.calc_component_stats:
                ao.print_components(
                    self.id,
                    self.time,
                    self.run_seed,
                    self.pop.pop_seed,
                    self.pop.connected_components(),
                    network_outdir,
                    self.params.classes.races,
                )

            if self.params.outputs.network.calc_network_stats:
                self.network_utils.write_network_stats(
                    network_outdir, self.id, self.time
                )

            if self.params.outputs.network.edge_list:
                self.network_utils.write_graph_edgelist(
                    network_outdir, self.id, self.time
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
            self.pop.all_agents, self.deaths, self.params, self.features, self.time,
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
                self.pop.hiv_agents.num_members(),
                sum([1 for a in self.pop.all_agents if a.incar.active]),  # type: ignore[attr-defined]
                sum([1 for a in self.pop.all_agents if a.high_risk.active]),  # type: ignore[attr-defined]
                sum([1 for a in self.pop.all_agents if a.prep.active]),  # type: ignore[attr-defined]
            )
        )

        self.timeline_scaling()

        self.update_all_agents()

        stats = ao.get_stats(
            self.pop.all_agents, self.deaths, self.params, self.features, self.time,
        )
        self.print_stats(stats, outdir)

        print(("Number of relationships: {}".format(len(self.pop.relationships))))
        self.pop.all_agents.print_subsets()

    def update_all_agents(self):
        """
        The core of the model.  For a time step, update all of the agents and relationships:

        1. Create an agent zero (if enabled and the time is right)
        2. Update partner assignments (create new relationships as needed)
        3. Agents in relationships interact
        4. Update features at the population level
        5. Update each agent's status for:
            * age
            * hiv
            * all features (agent level)
        6. End relationships with no remaining duration
        7. Agent death/replacement
        """
        # If agent zero enabled, create agent zero at the beginning of main loop.
        if (
            self.time == self.params.agent_zero.start_time
            and self.params.features.agent_zero
        ):
            self.make_agent_zero()

        if not self.params.features.static_network:
            self.pop.update_partner_assignments(t=self.time)
            if self.pop.enable_graph:
                self.pop.trim_graph()

        for rel in self.pop.relationships:
            self.agents_interact(rel)

        # TODO add check for whether feature is on somehow
        for feature in self.features:
            feature.update_pop(self)

        for agent in self.pop.all_agents:
            # happy birthday agents!
            if (
                self.time > 0
                and (self.time % self.params.model.time.steps_per_year) == 0
            ):
                agent.age += 1

            if agent.hiv:
                # If HIV hasn't started, ignore
                if self.time >= self.params.hiv.start_time:
                    self.diagnose_hiv(agent)
                    self.progress_to_aids(agent)

            for feature in self.features:
                agent_feature = getattr(agent, feature.name)
                agent_feature.update_agent(self)

        # If static network, ignore relationship progression
        if not self.params.features.static_network:
            for rel in copy(self.pop.relationships):
                if rel.progress():
                    self.pop.remove_relationship(rel)

        if self.params.features.die_and_replace:
            self.die_and_replace()

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
            self.hiv_convert(agent_zero)
        elif self.params.agent_zero.fallback and max_agent is not None:
            self.hiv_convert(max_agent)
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

    def agents_interact(self, rel: "ag.Relationship") -> bool:
        """
        Let an agent interact with a partner.

        Based on the interaction types of the relationship, interact in the following ways:

        * Peer Change Agent
        * Injection
        * Sex

        args:
            rel : The relationship that the agents interact in

        returns:
            whether the agents interacted
        """
        interaction_types = self.params.classes.bond_types[rel.bond_type].acts_allowed
        # If either agent is incarcerated, skip their interaction
        if rel.agent1.incar.active or rel.agent2.incar.active:  # type: ignore[attr-defined]
            return False

        agents_interacted = False
        for interaction_type in interaction_types:
            interaction = self.interactions[interaction_type]
            agents_interacted = interaction.interact(self, rel) or agents_interacted

        return agents_interacted

    def get_transmission_probability(self, interaction: str, agent, partner) -> float:
        """
        Determines the probability of a transmission event based on
            interaction type. For sex acts, transmission probability is a
            function of the acquisition probability of the HIV- agent's sex role
            and the HIV+ agent's haart adherence, acute status, and dx risk reduction

        args:
            interaction : "injection" or "sex"
            agent: HIV+ Agent
            partner: HIV- Agent

        returns:
            probability of transmission from agent to partner
        """
        # Logic for if needle or sex type interaction
        p: float
        assert interaction in ("injection", "sex"), (
            f"Invalid interaction type {interaction}. Only sex and injection acts "
            f"supported. "
        )

        # get baseline probabilities
        if interaction == "injection":
            p = self.params.partnership.injection.transmission.base
        elif interaction == "sex":
            agent_sex_role = agent.sex_role
            partner_sex_role = partner.sex_role

            # get partner's sex role during acts
            if partner_sex_role == "versatile":  # versatile partner takes
                # "opposite" position of agent
                if agent_sex_role == "insertive":
                    partner_sex_role = "receptive"
                elif agent_sex_role == "receptive":
                    partner_sex_role = "insertive"
                else:
                    partner_sex_role = "versatile"  # if both versatile, can switch
                    # between receptive and insertive by act

            # get probability of sex acquisition given HIV- partner's position
            p = partner.location.params.partnership.sex.acquisition[partner.sex_type][
                partner_sex_role
            ]

        # feature specific risk adjustment
        for feature in self.features:
            agent_feature = getattr(agent, feature.name)
            p *= agent_feature.get_transmission_risk_multiplier(self.time, interaction)

            partner_feature = getattr(partner, feature.name)
            p *= partner_feature.get_acquisition_risk_multiplier(self.time, interaction)

        # Scaling parameter for acute HIV infections
        if agent.get_acute_status(self.time):
            p *= agent.location.params.hiv.acute.infectivity

        # Scaling parameter for positively identified HIV agents
        if agent.hiv_dx:
            p *= 1 - agent.location.params.hiv.dx.risk_reduction[interaction]

        # Racial calibration parameter to attain proper race incidence disparity
        p *= partner.location.params.demographics[partner.race].hiv.transmission

        # Scaling parameter for per act transmission.
        p *= self.calibration.acquisition

        return p

    def hiv_convert(self, agent: "ag.Agent"):
        """
        Agent becomes HIV agent. Update all appropriate list and dictionaries.

        args:
            agent: The agent being converted
        """
        if not agent.hiv:
            agent.hiv = True
            agent.hiv_time = self.time
            agent.vaccine.active = False  # type: ignore[attr-defined]
            self.pop.hiv_agents.add_agent(agent)

        if agent.prep.active:  # type: ignore[attr-defined]
            agent.prep.progress(self, force=True)  # type: ignore[attr-defined]

    def diagnose_hiv(self, agent: "ag.Agent"):
        """
        Stochastically test the agent for HIV. If tested, mark the agent as diagnosed and trace their partners (if partner tracing enabled).

        args:
            agent: HIV positive agent to diagnose
        """
        sex_type = agent.sex_type
        race_type = agent.race
        diagnosed = agent.hiv_dx
        partner_tracing = agent.location.params.partner_tracing

        def diagnose(agent,):
            # agent's location's params used throughout as that is the agent who
            # would be interacting with the service
            agent.hiv_dx = True
            agent.hiv_dx_time = self.time
            self.pop.dx_counts[agent.race][agent.sex_type] += 1
            if (
                self.params.features.partner_tracing
                and partner_tracing.start_time <= self.time < partner_tracing.stop_time
            ):
                # Determine if each partner is found via partner tracing
                for ptnr in agent.get_partners(partner_tracing.bond_type):
                    if (
                        not ptnr.hiv_dx
                        and self.run_random.random() < partner_tracing.prob
                    ):
                        ptnr.partner_traced = True
                        ptnr.trace_time = self.time

        if not diagnosed:
            test_prob = agent.location.params.demographics[race_type][
                sex_type
            ].hiv.dx.prob

            # Rescale based on calibration param
            test_prob *= self.calibration.test_frequency

            if self.run_random.random() < test_prob:
                diagnose(agent)
            elif (
                agent.partner_traced
                and self.run_random.random() < partner_tracing.hiv.dx
                and self.time > agent.trace_time
            ):
                diagnose(agent)
        if self.time >= agent.trace_time + partner_tracing.trace_duration:
            # agents can only be traced during a specified period after their partner is
            # diagnosed. If past this time, remove ability to trace.
            agent.partner_traced = False

    def progress_to_aids(self, agent: "ag.Agent"):
        """
        Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        assert agent.hiv

        p = prob.adherent_prob(agent.haart.adherent) if agent.haart.active else 1  # type: ignore[attr-defined]

        if self.run_random.random() < p * agent.location.params.hiv.aids.prob:
            agent.aids = True

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
                    agent.hiv,
                    agent.aids,
                    agent.drug_type,
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
            # Remove agent from agent class and sub-sets
            self.pop.remove_agent(agent)

            new_agent = self.pop.create_agent(
                agent.location, agent.race, self.time, agent.sex_type, agent.drug_type
            )
            self.pop.add_agent(new_agent)
