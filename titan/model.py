# Imports
import random
from typing import Dict, List, Optional
from copy import copy
import os

import numpy as np  # type: ignore
import networkx as nx  # type: ignore
import nanoid  # type: ignore


from .agent import AgentSet, Agent, Relationship
from .population import Population
from .network import NetworkGraphUtils
from . import output as ao
from . import probabilities as prob
from . import utils
from .parse_params import ObjMap


class HIVModel:
    def __repr__(self):
        res = "\n"
        res += f"Seed: {self.run_seed}\n"
        res += f"Npop: {self.params.model.num_pop}\n"
        res += f"Time: {self.params.model.time.num_steps}\n"

        return res

    def __init__(
        self,
        params: ObjMap,
        population: Optional[Population] = None,
    ):
        """
        This is the core class used to simulate
            the spread of HIV and drug use in one geography.

        args:
            params: the parameter object for this model
            population: an initialized population to run the model on
        """

        self.params = params
        # pre-fetch commonly used param sub-sets for performance
        self.features = params.features
        self.calibration = params.calibration

        print("=== Begin Initialization Protocol ===\n")

        if population is None:
            print("\tGenerating new population")
            self.pop = Population(params)
        else:
            print("\tUsing provided population")
            self.pop = population

        self.network_utils: Optional[NetworkGraphUtils]
        if params.model.network.enable:
            self.network_utils = NetworkGraphUtils(self.pop.graph)
        else:
            self.network_utils = None

        print("\n\tCreating lists")
        # Other lists / dictionaries
        self.new_infections = AgentSet("new_infections")
        self.new_dx = AgentSet("new_dx")
        self.new_incar_release = AgentSet("new_incar_release")
        self.new_high_risk = AgentSet("new_high_risk")
        self.new_prep = AgentSet("new_prep")

        self.ssp_enrolled_risk = 0.0

        self.time = -1 * self.params.model.time.burn_steps  # burn is negative time
        self.id = nanoid.generate(size=8)

        # Set seed format. 0: pure random, else: fixed value
        self.run_seed = utils.get_check_rand_int(params.model.seed.run)
        print(f"\tRun seed was set to: {self.run_seed}")
        self.run_random = random.Random(self.run_seed)
        self.np_random = np.random.RandomState(self.run_seed)
        random.seed(self.run_seed)
        print(("\tFIRST RANDOM CALL {}".format(random.randint(0, 100))))

        print("\tResetting death count")
        self.deaths: List[Agent] = []  # Number of death

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
        self.new_infections.clear_set()
        self.new_dx.clear_set()
        self.new_high_risk.clear_set()
        self.new_incar_release.clear_set()
        self.new_prep.clear_set()
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
        if self.params.model.time.burn_steps > 0:
            print("\t===! Start Burn Loop !===")
        else:
            # make sure t0 things get printed
            stats = ao.get_stats(
                self.pop.all_agents,
                self.new_prep,
                self.new_infections,
                self.new_dx,
                self.new_high_risk,
                self.new_incar_release,
                self.deaths,
                self.params,
            )
            self.print_stats(stats, outdir)
        # burn is negative time, model run starts at t = 1
        for i in range(
            -1 * self.params.model.time.burn_steps, self.params.model.time.num_steps
        ):
            self.time += 1
            burn = True if self.time < 1 else False
            self.step(outdir, burn=burn)
            self.reset_trackers()

            if self.time == 0:
                if self.params.model.time.burn_steps > 0:
                    print("\t===! Burn Loop Complete !===")
                print("\t===! Start Main Loop !===")

        print("\t===! Main Loop Complete !===")

    def step(self, outdir: str, burn: bool = False):
        """
        A single time step in the model:

        1. Perform timeline_scaling updates to params if needed
        2. Update all agents
        3. Write/update reports with this timestep's data

        args:
            outdir: path to directory where reports should be saved
            burn: whether the model is in burn-in model (negative time)
        """
        print(f"\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t.: TIME {self.time}")
        print(
            "\tSTARTING HIV count:{}\tTotal Incarcerated:{}\tHR+:{}\t"
            "PrEP:{}".format(
                self.pop.hiv_agents.num_members(),
                sum([1 for a in self.pop.all_agents if a.incar]),
                self.pop.high_risk_agents.num_members(),
                sum([1 for a in self.pop.all_agents if a.prep]),
            )
        )

        self.timeline_scaling()

        self.update_all_agents(burn=burn)

        stats = ao.get_stats(
            self.pop.all_agents,
            self.new_prep,
            self.new_infections,
            self.new_dx,
            self.new_high_risk,
            self.new_incar_release,
            self.deaths,
            self.params,
        )
        self.print_stats(stats, outdir)

        print(("Number of relationships: {}".format(len(self.pop.relationships))))
        self.pop.all_agents.print_subsets()

    def update_all_agents(self, burn: bool = False):
        """
        The core of the model.  For a time step, update all of the agents and relationships:

        1. Create an agent zero (if enabled and the time is right)
        2. Update partner assignments (create new relationships as needed)
        3. Agents in relationships interact
        4. Update syringe services (if enabled)
        5. Update each agent's status for:
            * age
            * high risk
            * prep
            * incarceration
            * hiv
        6. End relationships with no remaining duration
        7. Agent death/replacement

        args:
            burn: whether the model is in burn-in period (negative time)
        """
        # If agent zero enabled, create agent zero at the beginning of main loop.
        if self.time == self.params.agent_zero.start_time and self.features.agent_zero:
            self.make_agent_zero()

        if not self.features.static_network:
            self.pop.update_partner_assignments(t=self.time)
            if self.pop.enable_graph:
                self.pop.trim_graph()

        for rel in self.pop.relationships:
            # If before hiv start time, ignore interactions
            if (self.time >= self.params.hiv.start) or (not burn and self.features.pca):
                self.agents_interact(rel)

        if self.features.syringe_services:
            self.update_syringe_services()

        for agent in self.pop.all_agents:
            # happy birthday agents!
            if (
                self.time > 0
                and (self.time % self.params.model.time.steps_per_year) == 0
            ):
                agent.age += 1

            if self.features.high_risk:
                self.update_high_risk(agent)

            if (
                self.features.pca
                and self.run_random.random()
                < agent.location.params.prep.pca.awareness.prob
                and not burn
            ):
                agent.prep_awareness = True
                if self.run_random.random() < agent.location.params.prep.pca.prep.prob:
                    self.initiate_prep(agent, force=True)

            if self.features.incar:
                self.incarcerate(agent)

            if (
                agent.msmw
                and self.run_random.random() < agent.location.params.msmw.hiv.prob
            ):
                self.hiv_convert(agent)

            if agent.hiv:
                agent.hiv_time += 1
                # If HIV hasn't started, ignore
                if self.time >= self.params.hiv.start:
                    self.diagnose_hiv(agent)
                    self.progress_to_aids(agent)

                    if self.features.haart:
                        self.update_haart(agent)
            else:
                if self.features.prep:
                    if self.time >= agent.location.params.prep.start:
                        if agent.prep:
                            self.discontinue_prep(agent)
                        elif agent.prep_eligible():
                            self.initiate_prep(agent)

                    if self.features.vaccine and not agent.prep:
                        self.advance_vaccine(
                            agent, vaxType=agent.location.params.vaccine.type, burn=burn
                        )

        if (
            self.features.prep
            and self.time == agent.location.params.prep.start
            and "RandomTrial" in agent.location.params.prep.target_model
        ):
            self.initialize_random_trial()

        # If static network, ignore relationship progression
        if not self.features.static_network:
            for rel in copy(self.pop.relationships):
                if rel.progress():
                    self.pop.remove_relationship(rel)

        if self.features.die_and_replace:
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
        if not self.features.timeline_scaling:
            return None

        # gather all of the param objects to be scaled
        params_set = {self.params}
        for location in self.pop.geography.locations.values():
            params_set.add(location.params)

        # iterate over each param and update the values if the time is right
        for params in params_set:
            for defn in params.timeline_scaling.timeline.values():
                param = defn.parameter
                if param != "ts_default":
                    if defn.time_start == self.time:
                        print(f"timeline scaling - {param}")
                        utils.scale_param(params, param, defn.scalar)
                    elif defn.time_stop == self.time:
                        print(f"timeline un-scaling - {param}")
                        utils.scale_param(params, param, 1 / defn.scalar)

    def update_high_risk(self, agent: Agent):
        """
        Update high risk agents or remove them from high risk pool
        """
        if agent not in self.pop.high_risk_agents:
            return None

        if agent.high_risk_time > 0:
            agent.high_risk_time -= 1
        else:
            self.pop.high_risk_agents.remove_agent(agent)
            agent.high_risk = False

            if self.features.incar:
                for bond in agent.location.params.high_risk.partnership_types:
                    agent.mean_num_partners[
                        bond
                    ] -= agent.location.params.high_risk.partner_scale
                    agent.mean_num_partners[bond] = max(
                        0, agent.mean_num_partners[bond]
                    )  # make sure not negative
                    agent.target_partners[bond] = utils.poisson(
                        agent.mean_num_partners[bond], self.np_random
                    )
                    while len(agent.partners[bond]) > agent.target_partners[bond]:
                        rel = utils.safe_random_choice(
                            agent.relationships, self.run_random
                        )
                        if rel is not None:
                            rel.progress(force=True)
                            self.pop.remove_relationship(rel)

    def initialize_random_trial(self):
        """
        Initialize a random trial in the population
        """
        assert (
            self.params.model.network.enable
        ), "Network must be enabled for random trial"

        print("Starting random trial")
        components = self.pop.connected_components()

        total_nodes = 0
        print(
            "Number of components",
            len([1 for comp in components if comp.number_of_nodes()]),
        )
        for comp in components:
            total_nodes += comp.number_of_nodes()
            if (
                self.run_random.random()
                < self.params.prep.random_trial.intervention.prob
            ):
                # Component selected as treatment pod!
                if not self.features.pca:
                    for ag in comp.nodes():
                        ag.random_trial_enrolled = True
                        if not ag.hiv and not ag.prep:
                            ag.intervention_ever = True
                            if (
                                self.run_random.random()
                                < ag.location.params.prep.target
                                and not ag.vaccine
                            ):
                                self.initiate_prep(ag, force=True)
                elif self.params.prep.pca.choice == "eigenvector":
                    centrality = nx.algorithms.centrality.eigenvector_centrality(comp)
                    assert len(centrality) >= 1, "Empty centrality"
                    ordered_centrality = sorted(centrality, key=centrality.get)
                    intervention_agent = False
                    for ag in ordered_centrality:
                        ag.random_trial_enrolled = True
                        if not ag.hiv and not intervention_agent:
                            ag.prep_awareness = True
                            ag.pca = True
                            ag.pca_suitable = True
                            ag.intervention_ever = True
                            intervention_agent = True

                    if not intervention_agent:
                        ag = ordered_centrality[0]
                elif self.params.prep.pca.choice == "bridge":
                    # list all edges that are bridges
                    for ag in comp.nodes:
                        ag.random_trial_enrolled = True

                    all_bridges = list(nx.bridges(comp))
                    comp_agents = [
                        agent
                        for agents in all_bridges
                        for agent in agents
                        if not agent.hiv
                    ]  # all suitable agents in bridges

                    if comp_agents:
                        chosen_agent = utils.safe_random_choice(
                            comp_agents, self.run_random
                        )  # select change agent
                        chosen_agent.prep_awareness = True  # make aware
                        chosen_agent.pca = True
                        chosen_agent.pca_suitable = True
                    else:
                        chosen_agent = list(comp.nodes)[0]
                        chosen_agent.pca = True

                elif self.params.prep.pca.choice == "random":
                    suitable_agent_choices = []
                    for ag in comp.nodes:
                        ag.random_trial_enrolled = True
                        if not ag.hiv:
                            suitable_agent_choices.append(ag)

                    if (
                        suitable_agent_choices
                    ):  # if there are agents who meet eligibility criteria,
                        # select one randomly
                        chosen_agent = utils.safe_random_choice(
                            suitable_agent_choices, self.run_random
                        )
                        chosen_agent.pca = True
                        chosen_agent.pca_suitable = True
                        chosen_agent.prep_awareness = True  # make aware
                        chosen_agent.intervention_ever = True
                    else:  # if no suitable agents, mark a non-suitable agent
                        chosen_agent = utils.safe_random_choice(
                            list(comp.nodes), self.run_random
                        )
                        chosen_agent.pca = True

        print(("Total agents in trial: ", total_nodes))

    def agents_interact(self, rel: Relationship) -> bool:
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
        if rel.agent1.incar or rel.agent2.incar:
            return False

        # Agent 1 is HIV+, Agent 2 is not, Agent 2 is succept
        if rel.agent1.hiv and not rel.agent2.hiv:
            agent = rel.agent1
            partner = rel.agent2
        # If Agent 2 is HIV and Agent 1 is not, Agent 1 is succept
        elif not rel.agent1.hiv and rel.agent2.hiv:
            agent = rel.agent2
            partner = rel.agent1
        else:  # neither agent is HIV or both are
            return False

        if "pca" in interaction_types and rel.duration < rel.total_duration:
            self.pca_interaction(rel)

        if "injection" in interaction_types:
            self.injection_transmission(agent, partner)

        if "sex" in interaction_types:
            self.sex_transmission(rel)

        return True

    def pca_interaction(self, rel: Relationship, force=False):
        """
        Simulate peer change agent interactions. Knowledge if one agent is aware and one unaware,
            opinion if one agent swaying the other.

        args:
            rel: The relationship PCA is happening in
            force: Whether to force knowledge dissemination and influence
        """

        assert (
            self.params.model.network.enable
        ), "Network must be enabled for pca interactions"

        def influence(agent, partner):
            agent_init_opinion = agent.prep_opinion
            partner_init_opinion = partner.prep_opinion
            agent_influence = nx.closeness_centrality(self.pop.graph, agent)
            partner_influence = nx.closeness_centrality(self.pop.graph, partner)

            if agent_influence > partner_influence:
                partner.prep_opinion = np.mean(
                    [agent.prep_opinion, partner.prep_opinion]
                )
            elif agent_influence < partner_influence:
                agent.prep_opinion = np.mean([agent.prep_opinion, partner.prep_opinion])

            if (
                agent_init_opinion
                < agent.location.params.prep.pca.opinion.threshold
                < agent.prep_opinion
            ):
                if self.run_random.random() < agent.location.params.prep.pca.prep.prob:
                    self.initiate_prep(agent, force=True)

            elif (
                partner_init_opinion
                < partner.location.params.prep.pca.opinion.threshold
                < partner.prep_opinion
            ):
                if (
                    self.run_random.random()
                    < partner.location.params.prep.pca.prep.prob
                ):
                    self.initiate_prep(partner, force=True)

        def knowledge_dissemination(partner):
            partner.prep_awareness = True
            if (
                partner.prep_opinion
                > partner.location.params.prep.pca.opinion.threshold
                and self.run_random.random()
                < partner.location.params.prep.pca.prep.prob
            ):
                self.initiate_prep(partner, force=True)

        def knowledge_transmission_probability():
            if rel.agent1.prep_awareness and rel.agent2.prep_awareness:
                p = self.params.prep.pca.opinion.transmission
            else:
                p = self.params.prep.pca.knowledge.transmission

            if num_acts == 1:
                p_total_transmission = p
            elif num_acts >= 1:
                p_total_transmission = 1.0 - utils.binom_0(num_acts, p)
            else:
                p_total_transmission = 0

            return p_total_transmission

        acts_prob = self.run_random.random()
        acts_bin = 0
        current_p_value = 0.0

        while acts_prob > current_p_value:
            acts_bin += 1
            current_p_value += self.params.partnership.pca.frequency[rel.bond_type][
                acts_bin
            ].prob
        min = self.params.partnership.pca.frequency[rel.bond_type][acts_bin].min
        max = self.params.partnership.pca.frequency[rel.bond_type][acts_bin].max
        if min == max:
            num_acts = min
        else:
            num_acts = self.run_random.randrange(min, max)

        if num_acts < 1 and not force:
            return
        if rel.agent1.prep_awareness and not rel.agent2.prep_awareness:
            if self.run_random.random() < knowledge_transmission_probability() or force:
                knowledge_dissemination(rel.agent2)
        elif not rel.agent1.prep_awareness and rel.agent2.prep_awareness:
            if self.run_random.random() < knowledge_transmission_probability() or force:
                knowledge_dissemination(rel.agent1)
        elif rel.agent1.prep_awareness and rel.agent2.prep_awareness or force:
            if self.run_random.random() < knowledge_transmission_probability() or force:
                influence(rel.agent1, rel.agent2)

    def injection_transmission(self, agent: Agent, partner: Agent):
        """
        Simulate random transmission of HIV between two PWID agents through injection.

        args:
            agent: PWID agent with HIV
            partner: PWID agent without HIV
        """

        assert agent.hiv
        assert not partner.hiv
        assert agent.drug_type == "Inj"
        assert partner.drug_type == "Inj"

        agent_race = agent.race
        agent_sex_type = agent.sex_type

        mean_num_acts = (
            agent.location.params.demographics[agent_race][
                agent_sex_type
            ].injection.num_acts
            * self.calibration.injection.act
        )
        share_acts = utils.poisson(mean_num_acts, self.np_random)

        if agent.ssp or partner.ssp:  # syringe services program risk
            p_unsafe_injection = self.ssp_enrolled_risk
        else:
            # If sharing, minimum of 1 share act
            if share_acts < 1:
                share_acts = 1

            p_unsafe_injection = agent.location.params.demographics[agent_race][
                agent_sex_type
            ].injection.unsafe_prob

            if agent.hiv_dx or partner.hiv_dx:  # diagnosis risk reduction
                p_unsafe_injection *= 1 - self.params.hiv.dx.risk_reduction.injection

        for n in range(share_acts):
            if self.run_random.random() > p_unsafe_injection:
                share_acts -= 1

        if share_acts >= 1.0:
            p = self.get_transmission_probability("injection", agent, partner)

            p_total_transmission: float
            if share_acts == 1:
                p_total_transmission = p
            else:
                p_total_transmission = 1.0 - utils.binom_0(share_acts, p)

            if self.run_random.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                self.hiv_convert(partner)

    def sex_transmission(self, rel: Relationship):
        """
        Simulate random transmission of HIV between two agents through Sex. One of the agents must have HIV.

        args:
            rel : Relationship
        """

        if rel.agent1.hiv:
            agent = rel.agent1
            partner = rel.agent2
        elif rel.agent2.hiv:
            agent = rel.agent2
            partner = rel.agent1
        else:
            raise ValueError("rel must have an agent with HIV")

        # HIV status of agent and partner
        # Everything from here is only run if one of them is HIV+
        if partner.hiv:
            return

        # unprotected sex probabilities for primary partnerships
        mean_sex_acts = (
            agent.get_number_of_sex_acts(self.run_random) * self.calibration.sex.act
        )
        total_sex_acts = utils.poisson(mean_sex_acts, self.np_random)

        # Get condom usage
        p_safe_sex = agent.location.params.demographics[agent.race][
            agent.sex_type
        ].safe_sex
        # increase condom usage if diagnosed
        if agent.hiv_dx or partner.hiv_dx:
            # Calculate probability of safe sex given risk reduction
            p_unsafe_sex = (1 - p_safe_sex) * (
                1 - self.params.hiv.dx.risk_reduction.sex
            )
            p_safe_sex *= 1 - p_unsafe_sex

        # Reduction of risk acts between partners for condom usage
        unsafe_sex_acts = total_sex_acts
        for n in range(unsafe_sex_acts):
            if self.run_random.random() < p_safe_sex:
                unsafe_sex_acts -= 1

        if unsafe_sex_acts >= 1:
            # agent is HIV+
            rel.total_sex_acts += unsafe_sex_acts
            p_per_act = self.get_transmission_probability("sex", agent, partner)

            p_total_transmission: float
            if unsafe_sex_acts == 1:
                p_total_transmission = p_per_act
            else:
                p_total_transmission = 1.0 - utils.binom_0(unsafe_sex_acts, p_per_act)

            if self.run_random.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                self.hiv_convert(partner)

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

        agent_sex_role = agent.sex_role
        partner_sex_role = partner.sex_role

        if interaction == "injection":
            p = self.params.partnership.injection.transmission.base
            if agent.haart:
                p *= agent.location.params.partnership.injection.transmission.haart_scaling[
                    agent.haart_adherence
                ].scale
        elif interaction == "sex":
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

            # scale based on HIV+ agent's haart status/adherence
            if agent.haart:
                p *= agent.location.params.partnership.sex.haart_scaling[
                    agent.sex_type
                ][agent.haart_adherence].prob

        # Scale if partner on PrEP
        if partner.prep:
            if partner.prep_type == "Oral":
                if partner.prep_adherence == 1:
                    p *= 1.0 - partner.location.params.prep.efficacy.adherent
                else:
                    p *= 1.0 - partner.location.params.prep.efficacy.non_adherant
            elif partner.prep_type == "Inj" and partner.prep_adherence == 1:
                p *= -1.0 * np.exp(-5.528636721 * partner.prep_load)

        # Scale if partner vaccinated
        if partner.vaccine:
            vaccine_type = partner.location.params.vaccine.type
            vaccine_time_months = (
                partner.vaccine_time / self.params.model.time.steps_per_year
            ) * 12

            if vaccine_type == "HVTN702":
                p *= np.exp(-2.88 + 0.76 * (np.log((vaccine_time_months + 0.001) * 30)))
            elif vaccine_type == "RV144":
                p *= np.exp(-2.40 + 0.76 * (np.log(vaccine_time_months)))

        # Scaling parameter for acute HIV infections
        if agent.get_acute_status(agent.location.params.hiv.acute.duration):
            p *= agent.location.params.hiv.acute.infectivity

        # Scaling parameter for positively identified HIV agents
        if agent.hiv_dx:
            p *= 1 - agent.location.params.hiv.dx.risk_reduction[interaction]

        # Tuning parameter for ART efficiency
        if agent.haart:
            p *= self.calibration.haart.transmission

        # Racial calibration parameter to attain proper race incidence disparity
        p *= partner.location.params.demographics[partner.race].hiv.transmission

        # Scaling parameter for per act transmission.
        p *= self.calibration.acquisition

        return p

    def hiv_convert(self, agent: Agent):
        """
        Agent becomes HIV agent. Update all appropriate list and dictionaries.

        args:
            agent: The agent being converted
        """
        if not agent.hiv:
            agent.hiv = True
            agent.hiv_time = 1
            agent.vaccine = False
            self.new_infections.add_agent(agent)
            self.pop.hiv_agents.add_agent(agent)

        if agent.prep:
            self.discontinue_prep(agent, force=True)

    def update_syringe_services(self):
        """
        Enroll PWID agents in syringe services
        """
        print(("\n\n!!!!Engaging syringe services program"))
        ssp_num_slots = 0
        ssp_agents = {agent for agent in self.pop.pwid_agents.members if agent.ssp}
        if self.features.syringe_services:
            for item in self.params.syringe_services.timeline.values():
                if item.time_start <= self.time < item.time_stop:
                    self.ssp_enrolled_risk = item.risk

                    ssp_num_slots = (item.num_slots_stop - item.num_slots_start) / (
                        item.time_stop - item.time_start
                    ) * (self.time - item.time_start) + item.num_slots_start

                    # If cap indicates all or no agents, do not change
                    # otherwise, find true number of slots through distribution
                    if 0 < ssp_num_slots < self.pop.pwid_agents.num_members():
                        ssp_num_slots = round(
                            self.run_random.betavariate(
                                ssp_num_slots,
                                self.pop.pwid_agents.num_members() - ssp_num_slots,
                            )
                            * self.pop.pwid_agents.num_members()
                        )
                    break

        target_set = utils.safe_shuffle(
            (self.pop.pwid_agents.members - ssp_agents), self.run_random
        )

        for agent in ssp_agents.copy():
            if len(ssp_agents) > ssp_num_slots:
                agent.ssp = False
                ssp_agents.remove(agent)

        if target_set:
            for agent in target_set:
                if len(ssp_agents) < ssp_num_slots:
                    agent.ssp = True
                    ssp_agents.add(agent)

        print(
            f"SSP has {ssp_num_slots} target slots with "
            f"{len(ssp_agents)} slots filled"
        )

    def become_high_risk(self, agent: Agent, duration: int = None):
        """
        Mark an agent as high risk and assign a duration to their high risk period

        args:
            agent: agent becoming high risk
            duration: duration of the high risk period, defaults to param value if not passed [params.high_risk.sex_based]
        """

        if not self.features.high_risk:
            return None

        if agent not in self.pop.high_risk_agents.members:
            self.pop.high_risk_agents.add_agent(agent)

        if not agent.high_risk_ever:
            self.new_high_risk.add_agent(agent)

        agent.high_risk = True
        agent.high_risk_ever = True

        if duration is not None:
            agent.high_risk_time = duration
        else:
            agent.high_risk_time = agent.location.params.high_risk.sex_based[
                agent.sex_type
            ].duration

    def incarcerate(self, agent: Agent):
        """
        Incarcerate an agent or update their incarceration variables

        args:
            agent: agent being updated
        """
        if not self.features.incar:
            return None

        hiv_bool = agent.hiv

        if hiv_bool:
            hiv_multiplier = agent.location.params.incar.hiv.multiplier
        else:
            hiv_multiplier = 1

        if agent.incar:
            agent.incar_time -= 1

            if agent.incar_time == 0:  # FREE AGENT
                self.new_incar_release.add_agent(agent)
                agent.incar = False
                if (
                    not agent.high_risk and self.features.high_risk
                ):  # If behavioral treatment on and agent HIV, ignore HR period.
                    self.become_high_risk(agent)
                    for bond in agent.location.params.high_risk.partnership_types:
                        agent.mean_num_partners[
                            bond
                        ] += agent.location.params.high_risk.partner_scale
                        agent.target_partners[bond] = utils.poisson(
                            agent.mean_num_partners[bond], self.np_random
                        )
                    self.pop.update_partnerability(agent)

                if hiv_bool:
                    if agent.haart:
                        if (
                            self.run_random.random()
                            <= agent.location.params.incar.haart.discontinue
                        ):  # 12% remain surpressed
                            agent.haart = False
                            agent.haart_adherence = 0

                        # END FORCE

        elif self.run_random.random() < (
            agent.location.params.demographics[agent.race][agent.sex_type].incar.prob
            * hiv_multiplier
            * self.calibration.incarceration
        ):
            incar_duration = agent.location.params.demographics[agent.race][
                agent.sex_type
            ].incar.duration.prob

            bin = current_p_value = 1
            p = self.run_random.random()
            while p >= current_p_value:
                current_p_value += incar_duration[bin].prob
                bin += 1

            timestay = self.run_random.randint(
                incar_duration[bin].min, incar_duration[bin].max
            )

            if hiv_bool:
                if not agent.hiv_dx:
                    if self.run_random.random() < agent.location.params.incar.hiv.dx:
                        agent.hiv_dx = True
                else:  # Then tested and HIV, check to enroll in ART
                    if (
                        self.run_random.random()
                        < agent.location.params.incar.haart.prob
                    ):
                        tmp_rnd = self.run_random.random()
                        haart_adh = agent.location.params.incar.haart.adherence
                        if tmp_rnd < haart_adh:
                            adherence = 5
                        else:
                            adherence = self.run_random.randint(1, 4)

                        # Add agent to HAART class set, update agent params
                        agent.haart = True
                        agent.haart_adherence = adherence
                        agent.haart_time = self.time

            agent.incar = True
            agent.incar_time = timestay

            # PUT PARTNERS IN HIGH RISK
            for bond in agent.location.params.high_risk.partnership_types:
                for partner in agent.partners[bond]:
                    if not partner.high_risk and self.features.high_risk:
                        if (
                            self.run_random.random()
                            < partner.location.params.high_risk.prob
                        ):
                            self.become_high_risk(partner)

    def diagnose_hiv(self, agent: Agent):
        """
        Stochastically test the agent for HIV. If tested, mark the agent as diagnosed and trace their partners (if partner tracing enabled).

        args:
            agent: HIV positive agent to diagnose
        """
        sex_type = agent.sex_type
        race_type = agent.race
        diagnosed = agent.hiv_dx
        partner_tracing = agent.location.params.partner_tracing

        def diagnose(
            agent,
        ):
            # agent's location's params used throughout as that is the agent who
            # would be interacting with the service
            agent.hiv_dx = True
            agent.hiv_dx_time = 1
            self.pop.dx_counts[agent.race][agent.sex_type] += 1
            self.new_dx.add_agent(agent)
            if (
                self.features.partner_tracing
                and partner_tracing.start <= self.time < partner_tracing.stop
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
        else:
            agent.hiv_dx_time += 1
        if self.time >= agent.trace_time + partner_tracing.trace_duration:
            # agents can only be traced during a specified period after their partner is
            # diagnosed. If past this time, remove ability to trace.
            agent.partner_traced = False

    def update_haart(self, agent: Agent):
        """
        Account for HIV treatment through highly active antiretroviral therapy
            (HAART).
            HAART was implemented in 1996, hence, there is treatment only after 1996.
            HIV treatment assumes that the agent knows their HIV+ status (`dx` is True).

        args:
            agent: agent being updated
        """
        if not self.features.haart:
            return None

        def initiate(agent):
            haart_adh = agent.location.params.demographics[agent.race][
                agent.sex_type
            ].haart.adherence
            if self.run_random.random() < haart_adh:
                adherence = 5
            else:
                adherence = self.run_random.randint(1, 4)

            # Add agent to HAART class set, update agent params
            agent.haart = True
            agent.haart_adherence = adherence
            agent.haart_time = self.time
            self.pop.haart_counts[agent.race][agent.sex_type] += 1

        # Check valid input
        assert agent.hiv

        # Determine probability of HIV treatment
        if agent.hiv_dx:
            haart_params = agent.location.params.demographics[agent.race][
                agent.sex_type
            ].haart
            agent_dx_time = agent.hiv_dx_time

            # Go on HAART
            if not agent.haart:
                if agent.location.params.hiv.haart_cap:
                    # if HAART is based on cap instead of prob, determine number of
                    # HAART agents based on % of diagnosed agents
                    num_dx_agents = self.pop.dx_counts[agent.race][agent.sex_type]
                    num_haart_agents = self.pop.haart_counts[agent.race][agent.sex_type]

                    if num_haart_agents < (
                        agent.location.params.demographics[agent.race][
                            agent.sex_type
                        ].haart.prob
                        * num_dx_agents
                    ):
                        initiate(agent)
                else:
                    if not agent.haart_ever:
                        for defn in agent.location.params.haart.start_haart:
                            if defn.start_time == agent.hiv_dx_time:
                                haart_prob = defn.prob
                                break
                    else:
                        haart_prob = haart_params.reinit_prob

                    haart_prob *= self.calibration.haart.coverage
                    if self.run_random.random() < haart_prob:
                        initiate(agent)
                            
            # Go off HAART
            elif agent.haart and self.run_random.random() < haart_params.discontinue:
                agent.haart = False
                agent.haart_adherence = 0
                agent.haart_time = 0
                self.pop.haart_counts[agent.race][agent.sex_type] -= 1

    def discontinue_prep(self, agent: Agent, force: bool = False):
        """
        Update agent's PrEP status and discontinue stochastically or if `force` is True

        args:
            agent: agent being updated
            force: whether to force discontinuation of PrEP
        """
        # Agent must be on PrEP to discontinue PrEP
        assert agent.prep

        # If force flag set, auto kick off prep.
        if force:
            self.pop.prep_counts[agent.race] -= 1
            agent.prep = False
            agent.prep_load = 0.0
            agent.prep_last_dose = 0
            return None

        # else if agent is on PrEP, see if they should discontinue
        if (
            self.run_random.random()
            < agent.location.params.demographics[agent.race][
                agent.sex_type
            ].prep.discontinue
            and agent.prep_type == "Oral"
        ):
            self.pop.prep_counts[agent.race] -= 1
            agent.prep = False
            agent.prep_type = ""

        if agent.prep_type == "Inj":
            agent.update_prep_load()
            # agent timed out of prep
            if not agent.prep:
                self.pop.prep_counts[agent.race] -= 1

    def advance_vaccine(self, agent: Agent, vaxType: str, burn: bool):
        """
        Progress vaccine. Agents may receive injection or progress in time
            since injection.

        args:
            agent: agent being updated
            vaxType: type of vaccine
            burn: whether the model is in burn-in mode
        """
        if not self.features.vaccine:
            return None

        vaccine_params = agent.location.params.demographics[agent.race][
            agent.sex_type
        ].vaccine

        if agent.vaccine and not burn:
            agent.vaccine_time += 1
            if (
                agent.location.params.vaccine.booster
                and agent.vaccine_time == vaccine_params.booster.interval
                and self.run_random.random() < vaccine_params.booster.prob
            ):
                agent.vaccinate(vaxType)

        elif self.time == agent.location.params.vaccine.start:
            if agent.location.params.vaccine.init == burn:  # both true or both false
                if self.run_random.random() < vaccine_params.prob:
                    agent.vaccinate(vaxType)

    def initiate_prep(self, agent: Agent, force: bool = False):
        """
        Place agents onto PrEP treatment. PrEP treatment assumes that the agent knows their HIV status is negative.

        args:
            agent : agent being updated
            force : whether to force the agent to enroll instead of using the appropriate algorithm per the prep params
        """

        def enroll_prep(self, agent: Agent):
            agent.enroll_prep(self.run_random)

            self.new_prep.add_agent(agent)
            self.pop.prep_counts[agent.race] += 1

        # agent must exist
        assert agent is not None

        # Prep only valid for agents not on prep and are HIV negative
        if agent.prep or agent.hiv:
            return

        # Determine probability of HIV treatment
        if force:
            enroll_prep(self, agent)
        else:
            if "Racial" in agent.location.params.prep.target_model:
                num_prep_agents = self.pop.prep_counts[agent.race]
                all_hiv_agents = self.pop.hiv_agents.members
                all_race = {a for a in self.pop.all_agents if a.race == agent.race}

                hiv_agents = len(all_hiv_agents & all_race)
                target_prep = (
                    len(all_race) - hiv_agents
                ) * agent.location.params.demographics[agent.race][
                    agent.sex_type
                ].prep.coverage

            else:
                num_prep_agents = sum(self.pop.prep_counts.values())
                target_prep = int(
                    (
                        self.pop.all_agents.num_members()
                        - self.pop.hiv_agents.num_members()
                    )
                    * agent.location.params.prep.target
                )

            if (
                num_prep_agents < target_prep
                and self.time >= agent.location.params.prep.start
                and agent.prep_eligible()
            ):
                enroll_prep(self, agent)

    def progress_to_aids(self, agent: Agent):
        """
        Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        assert agent.hiv

        p = prob.adherence_prob(agent.haart_adherence) if agent.haart else 1

        if self.run_random.random() < p * agent.location.params.hiv.aids.prob:
            agent.aids = True

    def die_and_replace(self):

        """
        Let agents die and replace the dead agent with a new agent randomly.
        """
        # die stage
        for agent in self.pop.all_agents:

            # agent incarcerated, don't evaluate for death
            if agent.incar:
                continue

            # death rate per 1 person-month
            p = (
                prob.get_death_rate(
                    agent.hiv,
                    agent.aids,
                    agent.drug_type,
                    agent.haart_adherence,
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
                agent.location, agent.race, self.time, agent.sex_type
            )
            self.pop.add_agent(new_agent)
