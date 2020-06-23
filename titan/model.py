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
    """
    :Purpose:
        This is the core class used to simulate
        the spread of HIV and drug use in one MSA
        (Metropolitan Statistical Area).

    :Input:
        params: ObjMap - the parameter object for this model
    """

    def __repr__(self):
        res = "\n"
        res += f"Seed: {self.run_seed}\n"
        res += f"Npop: {self.params.model.num_pop}\n"
        res += f"Time: {self.params.model.time.num_steps}\n"

        return res

    def __init__(self, params: ObjMap, population: Optional[Population] = None):

        self.params = params
        # pre-fetch commonly used param sub-sets for performance
        self.features = params.features
        self.prep = params.prep
        self.demographics = params.demographics
        self.calibration = params.calibration
        self.high_risk = params.high_risk
        self.vaccine = params.vaccine
        self.incar = params.incar

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

    def print_stats(self, stat: Optional[Dict[str, Dict[str, int]]], outdir: str):
        if stat is not None:
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
                    network_outdir, curtime=self.time, label=f"{self.id}",
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
        Core of the model:
            1. Prints networkReport for first agents.
            2. Makes agents become HIV (used for current key_time tracking for acute)
            3. Loops over all time steps
                a. _update AllAgents()
                b. reset death count
                c. _ self.die_and_replace()
                d. self._update_population()
                e. self._reset_partner_count()
        """
        if self.params.model.time.burn_steps > 0:
            print("\t===! Start Burn Loop !===")
        else:
            # make sure t0 things get printed
            self.print_stats(None, outdir)
        # burn is negative time, model run starts at t = 1
        for i in range(
            -1 * self.params.model.time.burn_steps, self.params.model.time.num_steps
        ):
            self.time += 1
            burn = True if self.time < 0 else False
            self.step(outdir, burn=burn)
            self.reset_trackers()

            if self.time == 0:
                if self.params.model.time.burn_steps > 0:
                    print("\t===! Burn Loop Complete !===")
                print("\t===! Start Main Loop !===")

        print("\t===! Main Loop Complete !===")

    def step(self, outdir: str, burn: bool = False):
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
        :Purpose:
            Update agents.  For a time step, update all of the agents and relationships

        :Input:
            agent, time

        :Output:
            none
        """
        # If agent zero enabled, create agent zero at the beginning of main loop.
        if self.time == self.params.agent_zero.start_time and self.features.agent_zero:
            self.make_agent_zero()

        if not self.features.static_network:
            self.pop.update_partner_assignments(t=self.time)
            if self.pop.enable_graph:
                self.pop.initialize_graph()

        for rel in self.pop.relationships:
            # If in burn, ignore interactions
            if self.time > self.params.hiv.start:
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
                and self.run_random.random() < self.prep.pca.awareness.prob
                and not burn
            ):
                agent.prep_awareness = True
                if self.run_random.random() < self.prep.pca.prep.prob:
                    self.initiate_prep(agent, force=True)

            if self.features.incar:
                self.incarcerate(agent)

            if agent.msmw and self.run_random.random() < self.params.msmw.hiv.prob:
                self.hiv_convert(agent)

            if agent.hiv:
                # If HIV hasn't started, ignore
                if self.time > self.params.hiv.start:
                    self.diagnose_hiv(agent)
                    self.progress_to_aids(agent)

                    if self.features.haart:
                        self.update_haart(agent)
                        agent.hiv_time += 1
            else:
                if self.features.prep:
                    if self.time >= self.prep.start:
                        if agent.prep:
                            self.discontinue_prep(agent)
                        elif (
                            agent.prep_eligible(
                                self.prep.target_model,
                                self.params.partnership.ongoing_duration,
                            )
                            and self.prep.target_model != "RandomTrial"
                        ):
                            self.initiate_prep(agent)

                    if self.features.vaccine and not agent.prep:
                        self.advance_vaccine(
                            agent, vaxType=self.vaccine.type, burn=burn
                        )

        if (
            self.features.prep
            and self.time == self.prep.start
            and self.prep.target_model == "RandomTrial"
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
        bond_type = self.params.agent_zero.bond_type
        interaction_type = self.params.agent_zero.interaction_type
        bonds = [
            i
            for i in self.params.classes.bond_types.values()
            if interaction_type in i.acts_allowed
        ]
        print(bonds)
        zero_eligible = [
            agent
            for agent in self.pop.all_agents.members
            if len(agent.partners[bond_type]) >= self.params.agent_zero.num_partners
        ]
        agent_zero = utils.safe_random_choice(zero_eligible, self.run_random)
        if agent_zero:
            self.hiv_convert(agent_zero)
        else:
            raise ValueError("No agent zero!")

    def update_high_risk(self, agent: Agent):
        """
        :Purpose:
            Update high risk agents or remove them from high risk pool
        """
        if agent not in self.pop.high_risk_agents:
            return None

        if agent.high_risk_time > 0:
            agent.high_risk_time -= 1
            if (
                agent.so == "HM"
                and self.features.prep
                and (self.prep.target_model in ("high_risk", "incarcerated_high_risk"))
            ):
                for part in agent.iter_partners():
                    if not (part.hiv or part.vaccine):
                        self.initiate_prep(part)
        else:
            self.pop.high_risk_agents.remove_agent(agent)
            agent.high_risk = False

            if self.features.incar:
                for bond in self.params.high_risk.partnership_types:
                    agent.mean_num_partners[bond] -= self.high_risk.partner_scale
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
        :Purpose:
            Initialize random trial in population
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
                        ag.intervention_comp = True
                        if not ag.hiv and not ag.prep:
                            ag.intervention_ever = True
                            if (
                                self.run_random.random() < self.prep.target
                                and not ag.vaccine
                            ):
                                self.initiate_prep(ag, force=True)
                elif self.prep.pca.choice == "eigenvector":
                    centrality = nx.algorithms.centrality.eigenvector_centrality(comp)
                    assert len(centrality) >= 1, "Empty centrality"
                    ordered_centrality = sorted(centrality, key=centrality.get)
                    intervention_agent = False
                    for ag in ordered_centrality:
                        ag.intervention_comp = True
                        if not ag.hiv:
                            ag.prep_awareness = True
                            ag.pca = True
                            ag.pca_suitable = True
                            ag.intervention_ever = True
                            intervention_agent = True
                            break
                    if not intervention_agent:
                        ag = ordered_centrality[0]
                elif self.prep.pca.choice == "bridge":
                    # list all edges that are bridges
                    for ag in comp.nodes:
                        ag.intervention_comp = True

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

                elif self.prep.pca.choice == "random":
                    suitable_agent_choices = []
                    for ag in comp.nodes:
                        ag.intervention_comp = True
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

    def agents_interact(self, rel: Relationship):
        """
        :Purpose:
            Let PWID agent interact with a partner.
            Update PWID agents:
                1 - determine transition type
                2 - Injection rules
                3 - Sex rules
                4 - HIV transmission

        :Input:

            rel : Relationship

            rand_gen : random number generator

        Output:
            boolean : whether interaction happened

        """
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

        interaction_types = self.params.classes.bond_types[rel.bond_type].acts_allowed

        if "pca" in interaction_types and rel.duration < rel.total_duration:
            self.pca_interaction(rel)

        if "injection" in interaction_types:
            self.injection_transmission(agent, partner)

        if "sex" in interaction_types:
            self.sex_transmission(rel)

        return True

    def pca_interaction(self, rel: Relationship, force=False):
        """
        :Purpose:
            Simulate peer change agent interactions
            Knowledge if one agent is aware and one unaware,
            opinion if one agent swaying the other
        :Input:
            agent: Agent
            partner: Agent
            PCAtype: str, either 'Knowledge' or 'Opinion'
        :Output: -
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

            if self.run_random.random() < self.prep.pca.prep.prob:
                if (
                    agent_init_opinion
                    < self.prep.pca.opinion.threshold
                    < agent.prep_opinion
                ):
                    self.initiate_prep(agent, force=True)
                elif (
                    partner_init_opinion
                    < self.prep.pca.opinion.threshold
                    < partner.prep_opinion
                ):
                    self.initiate_prep(partner, force=True)

        def knowledge_dissemination(partner):
            partner.prep_awareness = True
            if (
                partner.prep_opinion > self.prep.pca.opinion.threshold
                and self.run_random.random() < self.prep.pca.prep.prob
            ):
                self.initiate_prep(partner, force=True)

        def knowledge_transmission_probability():
            if rel.agent1.prep_awareness and rel.agent2.prep_awareness:
                p = self.prep.pca.opinion.transmission
            else:
                p = self.prep.pca.knowledge.transmission

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
            current_p_value += self.params.partnership.interaction[rel.bond_type][
                acts_bin
            ].prob

        min = self.params.partnership.interaction[rel.bond_type][acts_bin].min
        max = self.params.partnership.interaction[rel.bond_type][acts_bin].max
        if min == max:
            num_acts = min
        else:
            num_acts = self.run_random.randrange(min, max)

        if num_acts < 1:
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
        :Purpose:
            Simulate random transmission of HIV between two PWID agents
            through injection.
            Agent must by HIV+ and partner not.

        :Input:
            agents : int
            partner : int
        :Output: -
        """

        assert agent.hiv
        assert not partner.hiv
        assert agent.drug_use == "Inj"
        assert partner.drug_use == "Inj"

        agent_race = agent.race
        agent_sex_type = agent.so

        mean_num_acts = (
            self.demographics[agent_race][agent_sex_type].injection.num_acts
            * self.calibration.injection.act
        )
        share_acts = utils.poisson(mean_num_acts, self.np_random)

        if agent.ssp:  # syringe services program risk
            p_unsafe_injection = self.ssp_enrolled_risk
        else:
            # If sharing, minimum of 1 share act
            if share_acts < 1:
                share_acts = 1

            p_unsafe_injection = self.demographics[agent_race][
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
        :Purpose:
            Simulate random transmission of HIV between two agents through Sex.
            Needed for all users. Sex is not possible in case the agent and
            assigned partner have incompatible Sex behavior. Given other logic,
            only one member of the relationship (the agent) has HIV at this time.

        :Input:
            rel : Relationship

        :Output:
            none
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
            agent.get_number_of_sex_acts(self.np_random, self.params)
            * self.calibration.sex.act
        )
        total_sex_acts = utils.poisson(mean_sex_acts, self.np_random)

        # Get condom usage
        p_safe_sex = self.demographics[agent.race][agent.so].safe_sex
        # increase condom usage if diagnosed
        if agent.hiv_dx or partner.hiv_dx:
            p_safe_sex *= 1 - self.params.hiv.dx.risk_reduction.sex

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
        """ Decriptor
        :Purpose:
            Determines the probability of a transmission event based on
            interaction type. For sex acts, transmission probability is a
            function of the acquisition probability of the HIV- agent's sex role
            and the HIV+ agent's haart adherence, acute status, and dx risk reduction

        :Input:
            interaction : str - "injection" or "sex"

        :Output:
            probability : float
        """
        # Logic for if needle or sex type interaction
        p: float
        assert interaction in ("injection", "sex",), (
            f"Invalid interaction type {interaction}. Only sex and injection acts "
            f"supported. "
        )

        agent_sex_role = agent.sex_role
        partner_sex_role = partner.sex_role

        if interaction == "injection":
            p = self.params.partnership.injection.transmission.base
            if agent.haart:
                p *= self.params.partnership.injection.transmission.haart_scaling[
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

            p = self.params.partnership.sex.acquisition[partner.so][partner_sex_role]

            # scale based on HIV+ agent's haart status/adherence
            if agent.haart:
                p *= self.params.partnership.sex.haart_scaling[agent.so][
                    agent.haart_adherence
                ].prob

        # Scale if partner on PrEP
        if partner.prep:
            if partner.prep_type == "Oral":
                if partner.prep_adherence == 1:
                    p *= 1.0 - self.prep.efficacy.adherent
                else:
                    p *= 1.0 - self.prep.efficacy.non_adherant
            elif partner.prep_type == "Inj" and partner.prep_adherence == 1:
                p *= -1.0 * np.exp(-5.528636721 * partner.prep_load)

        # Scale if partner vaccinated
        if partner.vaccine:
            assert self.vaccine.type in [
                "HVTN702",
                "RV144",
            ], f"Vaccine type {self.vaccine.type} not recognized"
            vaccine_time_months = (
                partner.vaccine_time / self.params.model.time.steps_per_year
            ) * 12
            if self.vaccine.type == "HVTN702":
                p *= np.exp(-2.88 + 0.76 * (np.log((vaccine_time_months + 0.001) * 30)))
            elif self.vaccine.type == "RV144":
                p *= np.exp(-2.40 + 0.76 * (np.log(vaccine_time_months)))

        # Scaling parameter for acute HIV infections
        if agent.get_acute_status(self.params.hiv.acute.duration):
            p *= self.params.hiv.acute.infectivity

        # Scaling parameter for positively identified HIV agents
        if agent.hiv_dx:
            p *= 1 - self.params.hiv.dx.risk_reduction[interaction]

        # Tuning parameter for ART efficiency
        if agent.haart:
            p *= self.params.calibration.haart.transmission

        # Racial calibration parameter to attain proper race incidence disparity
        p *= self.params.demographics[partner.race].hiv.transmission

        # Scaling parameter for per act transmission.
        p *= self.params.calibration.acquisition

        return p

    def hiv_convert(self, agent: Agent):
        """
        :Purpose:
            agent becomes HIV agent. Update all appropriate list and
            dictionaries.

        :Input:
            agent : int
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
        :Purpose:
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
            agent.high_risk_time = self.high_risk.sex_based[agent.so].duration

    def incarcerate(self, agent: Agent):
        """
        :Purpose:
            To incarcerate an agent or update their incarceration variables

        :Input:
            agent : int
        """
        if not self.features.incar:
            return None

        hiv_bool = agent.hiv

        if hiv_bool:
            hiv_multiplier = self.incar.hiv.multiplier
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
                    for bond in self.params.high_risk.partnership_types:
                        agent.mean_num_partners[bond] += self.high_risk.partner_scale
                        agent.target_partners[bond] = utils.poisson(
                            agent.mean_num_partners[bond], self.np_random
                        )
                    self.pop.update_partnerability(agent)

                if hiv_bool:
                    if agent.haart:
                        if (
                            self.run_random.random() <= self.incar.haart.discontinue
                        ):  # 12% remain surpressed
                            agent.haart = False
                            agent.haart_adherence = 0

                        # END FORCE

        elif self.run_random.random() < (
            self.demographics[agent.race][agent.so].incar.prob
            * hiv_multiplier
            * self.calibration.incarceration
        ):
            incar_duration = self.demographics[agent.race][agent.so].incar.duration.prob

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
                    if self.run_random.random() < self.incar.hiv.dx:
                        agent.hiv_dx = True
                else:  # Then tested and HIV, check to enroll in ART
                    if self.run_random.random() < self.incar.haart.prob:
                        tmp_rnd = self.run_random.random()
                        haart_adh = self.incar.haart.adherence
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
            for bond in self.params.high_risk.partnership_types:
                for partner in agent.partners[bond]:
                    if not partner.high_risk and self.features.high_risk:
                        if self.run_random.random() < self.high_risk.prob:
                            self.become_high_risk(partner)

                    if self.features.prep and (
                        self.prep.target_model in ("Incar", "IncarHR")
                    ):
                        # Attempt to put partner on prep if less than probability
                        if not partner.hiv and not agent.vaccine:
                            self.initiate_prep(partner)

    def diagnose_hiv(self, agent: Agent):
        """
        :Purpose:
            Test the agent for HIV. If detected, add to identified list.

        :Input:
            agent : agent_Class

        :Output:
            none
        """
        sex_type = agent.so
        race_type = agent.race
        diagnosed = agent.hiv_dx

        def diagnose(agent):
            agent.hiv_dx = True
            self.pop.num_dx_agents += 1
            self.new_dx.add_agent(agent)
            if (
                self.features.partner_tracing
                and self.params.partner_tracing.start
                <= self.time
                < self.params.partner_tracing.stop
            ):
                # Determine if each partner is found via partner tracing
                for bond in self.params.partner_tracing.bond_type:
                    for ptnr in agent.partners.get(bond, []):
                        if (
                            ptnr.hiv
                            and not ptnr.hiv_dx
                            and self.run_random.random()
                            < self.params.partner_tracing.prob
                        ):
                            ptnr.partner_traced = True
                            ptnr.trace_time = self.time

        if not diagnosed:
            test_prob = self.demographics[race_type][sex_type].hiv.dx.prob

            if (
                self.params.agent_zero.start_time
                <= self.time
                < self.params.agent_zero.start_time + self.params.agent_zero.dx_time
            ):
                test_prob *= self.params.agent_zero.dx_scalar

            # Rescale based on calibration param
            test_prob *= self.calibration.test_frequency

            if self.run_random.random() < test_prob:
                diagnose(agent)
            elif (
                agent.partner_traced
                and self.run_random.random() < self.params.partner_tracing.hiv.dx
                and self.time > agent.trace_time
            ):
                diagnose(agent)
        if self.time >= agent.trace_time + self.params.partner_tracing.trace_time:
            # agents can only be traced during a specified period after their partner is
            # diagnosed. If past this time, remove ability to trace.
            agent.partner_traced = False

    def update_haart(self, agent: Agent):
        """
        :Purpose:
            Account for HIV treatment through highly active antiretroviral therapy
            (HAART).
            HAART was implemented in 1996, hence, there is treatment only after 1996.
            HIV treatment assumes that the agent knows their HIV+ status.

        :Input:
            agent : Agent

        :Output:
            none
        """
        if not self.features.haart:
            return None

        def initiate(agent):
            haart_adh = self.demographics[agent_race][agent_so].haart.adherence
            if self.run_random.random() < haart_adh:
                adherence = 5
            else:
                adherence = self.run_random.randint(1, 4)

            # Add agent to HAART class set, update agent params
            agent.haart = True
            agent.haart_adherence = adherence
            agent.haart_time = self.time
            self.pop.num_haart_agents += 1

        # Check valid input
        assert agent.hiv

        agent_haart = agent.haart
        agent_race = agent.race
        agent_so = agent.so

        # Determine probability of HIV treatment
        if agent.hiv_dx:
            # Go on HAART
            if not agent_haart:
                if self.params.hiv.haart_cap:
                    # if HAART is based on cap instead of prob, determine number of
                    # HAART agents based on % of diagnosed agents
                    if (
                        self.pop.num_haart_agents
                        < self.demographics[agent_race][agent_so].haart.prob
                        * self.pop.num_dx_agents
                    ):
                        initiate(agent)
                else:
                    if self.run_random.random() < (
                        self.demographics[agent_race][agent_so].haart.prob
                        * self.calibration.haart.coverage
                    ):
                        initiate(agent)
            # Go off HAART
            elif (
                agent_haart
                and self.run_random.random()
                < self.demographics[agent_race][agent_so].haart.discontinue
            ):
                agent.haart = False
                agent.haart_adherence = 0
                agent.haart_time = 0
                self.pop.num_haart_agents -= 1

    def discontinue_prep(self, agent: Agent, force: bool = False):
        # Agent must be on PrEP to discontinue PrEP
        assert agent.prep

        # If force flag set, auto kick off prep.
        if force:
            self.pop.prep_counts[agent.race] -= 1
            agent.prep = False
            agent.prep_reason = []
            agent.prep_load = 0.0
            agent.prep_last_dose = 0
            return None

        # else if agent is on PrEP, see if they should discontinue
        if (
            self.run_random.random()
            < self.demographics[agent.race][agent.so].prep.discontinue
            and agent.prep_type == "Oral"
        ):
            self.pop.prep_counts[agent.race] -= 1
            agent.prep = False
            agent.prep_type = ""
            agent.prep_reason = []

        if agent.prep_type == "Inj":
            agent.update_prep_load(self.params)
            # agent timed out of prep
            if not agent.prep:
                self.pop.prep_counts[agent.race] -= 1

    def advance_vaccine(self, agent: Agent, vaxType: str, burn: bool):
        """
        :Purpose:
            Progress vaccine. Agents may receive injection or progress in time
            since injection.

        :Input:
            agent: Agent

        :Output:
            none
        """
        if not self.features.vaccine:
            return None

        if agent.vaccine and not burn:
            agent.vaccine_time += 1
            if (
                self.vaccine.booster
                and agent.vaccine_time
                == self.demographics[agent.race][agent.so].vaccine.booster.interval
                and self.run_random.random()
                < self.demographics[agent.race][agent.so].vaccine.booster.prob
            ):
                agent.vaccinate(vaxType)

        elif self.time == self.vaccine.start:
            if self.vaccine.init == burn:  # both true or both false
                if (
                    self.run_random.random()
                    < self.demographics[agent.race][agent.so].vaccine.prob
                ):
                    agent.vaccinate(vaxType)

    def initiate_prep(self, agent: Agent, force: bool = False):
        """
        :Purpose:
            Place agents onto PrEP treatment.
            PrEP treatment assumes that the agent knows their HIV+ status is negative.

        :Input:
            agent : Agent
            force : default is `False`

        :Output:
            none
        """

        def enroll_prep(self, agent: Agent):
            agent.enroll_prep(self.params, self.run_random)

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
            if self.prep.target_model == "Racial":
                num_prep_agents = self.pop.prep_counts[agent.race]
            else:
                num_prep_agents = sum(self.pop.prep_counts.values())

            if self.prep.target_model in ("Incar", "IncarHR"):
                if self.run_random.random() < self.prep.target:
                    enroll_prep(self, agent)
                return None
            elif self.prep.target_model == "Racial":
                all_hiv_agents = self.pop.hiv_agents.members
                all_race = {a for a in self.pop.all_agents if a.race == agent.race}

                hiv_agents = len(all_hiv_agents & all_race)
                target_prep = (len(all_race) - hiv_agents) * self.demographics[
                    agent.race
                ][agent.so].prep.coverage

            else:
                target_prep = int(
                    (
                        self.pop.all_agents.num_members()
                        - self.pop.hiv_agents.num_members()
                    )
                    * self.prep.target
                )

            if self.prep.target_model in ("Incar", "IncarHR"):
                if self.run_random.random() < self.prep.target:
                    enroll_prep(self, agent)
            elif (
                num_prep_agents < target_prep
                and self.time >= self.prep.start
                and agent.prep_eligible(
                    self.prep.target_model, self.params.partnership.ongoing_duration
                )
            ):
                enroll_prep(self, agent)

    def progress_to_aids(self, agent: Agent):
        """
        :Purpose:
            Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        assert agent.hiv

        if not agent.haart:
            p = prob.adherence_prob(agent.haart_adherence)

            if self.run_random.random() < p * self.params.hiv.aids.prob:
                agent.aids = True

    def die_and_replace(self):

        """
        :Purpose:
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
                    agent.drug_use,
                    agent.haart_adherence,
                    self.demographics[agent.race],
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

            new_agent = self.pop.create_agent(agent.race, agent.so)
            self.pop.add_agent(new_agent)
