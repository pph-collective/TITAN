#!/usr/bin/env python3
# encoding: utf-8

# Imports
import random
from typing import Dict, List, Sequence, Optional
import uuid

import numpy as np  # type: ignore
from scipy.stats import binom  # type: ignore
from scipy.stats import poisson  # type: ignore
import networkx as nx  # type: ignore
from dotmap import DotMap  # type: ignore


from .agent import AgentSet, Agent, Relationship
from .population import Population
from .network import NetworkGraphUtils
from . import output as ao
from . import probabilities as prob
from .partnering import sex_possible
from . import utils


class HIVModel:
    """
    :Purpose:
        This is the core class used to simulate
        the spread of HIV and drug use in one MSA
        (Metropolitan Statistical Area).

    :Input:
        params: DotMap - the parameter object for this model
    """

    def __repr__(self):
        res = "\n"
        res += "Seed: %d\n" % (self.run_seed)
        res += "Npop: %d\n" % (self.params.model.num_pop)
        res += "Time: %d\n" % (self.params.model.time_range)

        return res

    def __init__(self, params: DotMap, population: Optional[Population] = None):

        self.params = params

        self.run_seed = utils.get_check_rand_int(params.model.seed.run)

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

        self.total_dx = 0
        self.needle_exchange = False

        self.prep_agents: Dict[str, Dict[str, int]] = {}
        for race in params.classes.races:
            self.prep_agents[race] = {}
            for st in params.classes.sex_types:
                self.prep_agents[race][st] = 0

        # Set seed format. 0: pure random, -1: Stepwise from 1 to nRuns, else: fixed value
        print(("\tRun seed was set to:", self.run_seed))
        self.run_random = random.Random(self.run_seed)
        random.seed(self.run_seed)
        print(("\tFIRST RANDOM CALL %d" % random.randint(0, 100)))

        print("\tResetting death count")
        self.deaths: List[Agent] = []  # Number of death

        print("\n === Initialization Protocol Finished ===")

    def run(self, outdir):
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

        def print_stats(stat: Dict[str, Dict[str, int]], run_id: uuid.UUID):
            for report in self.params.outputs.reports:
                printer = getattr(ao, report)
                printer(
                    run_id,
                    t,
                    self.run_seed,
                    self.pop.pop_seed,
                    stat,
                    self.params,
                    outdir,
                )

        def reset_trackers():
            self.new_infections.clear_set()
            self.new_dx.clear_set()
            self.new_high_risk.clear_set()
            self.new_incar_release.clear_set()
            self.new_prep.clear_set()
            self.deaths = []

        def burn_simulation(duration: int):
            print(("\n === Burn Initiated for {} timesteps ===".format(duration + 1)))
            burn_conversions = 0
            for t in range(0, duration + 1):
                self.update_all_agents(t, burn=True)

                if self.params.features.die_and_replace:
                    self.die_and_replace()

                burn_conversions += self.new_infections.num_members()

                reset_trackers()

            self.pop.all_agents.print_subsets()

            print(f"\tBurn Cuml Inc:\t{burn_conversions}")

            print(" === Simulation Burn Complete ===")

        def make_agent_zero(num_partners: int):
            agent_zero = self.run_random.choice(self.pop.drug_use_inj_agents.members)
            for i in range(num_partners):
                self.pop.update_agent_partners(agent_zero)
            self.hiv_convert(agent_zero)

        run_id = uuid.uuid4()

        burn_simulation(self.params.model.burn_duration)

        print("\n === Begin Simulation Run ===")
        if self.params.outputs.network.draw_figures:
            self.network_utils.visualize_network(
                curtime=0, label="Seed" + str(self.run_seed),
            )

        if self.params.outputs.network.calc_component_stats:
            ao.print_components(
                run_id,
                0,
                self.run_seed,
                self.pop.pop_seed,
                self.pop.connected_components(),
                outdir,
            )

        print("\t===! Start Main Loop !===")

        # dictionary to hold results over time
        stats = {}

        # If we are using an agent zero method, create agent zero.
        if self.params.features.agent_zero:
            make_agent_zero(4)

        if self.params.outputs.network.edge_list:
            path = os.path.join(outdir, "network", "Edgelist_t0.txt")
            self.network_utils.write_graph_edgelist(path)

        for t in range(1, self.params.model.time_range + 1):
            print(f"\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t.: TIME {t}")
            if (
                self.params.outputs.network.draw_figures
                and t % self.params.outputs.print_frequency == 0
            ):
                self.network_utils.visualize_network(
                    curtime=t, label="Seed" + str(self.run_seed),
                )
            # todo: GET THIS TO THE NEW HIV COUNT

            print(
                "\tSTARTING HIV count:{}\tTotal Incarcerated:{}\tHR+:{}\tPrEP:{}".format(
                    self.pop.hiv_agents.num_members(),
                    self.pop.incarcerated_agents.num_members(),
                    self.pop.high_risk_agents.num_members(),
                    self.pop.intervention_prep_agents.num_members(),
                )
            )

            self.update_all_agents(t)

            if self.params.features.die_and_replace:
                self.die_and_replace()

            stats[t] = ao.get_stats(
                self.pop.all_agents,
                self.new_prep,
                self.new_infections,
                self.new_dx,
                self.pop.relationships,
                self.new_high_risk,
                self.new_incar_release,
                self.deaths,
                self.params,
            )
            print_stats(stats[t], run_id)

            print(("Number of relationships: %d" % len(self.pop.relationships)))
            self.pop.all_agents.print_subsets()

            self.total_dx += len(self.new_dx.members)
            if (
                self.total_dx > self.params.needle_exchange.init_at_pop
                and not self.needle_exchange
            ):
                self.enroll_needle_exchange()

            # RESET counters for the next time step
            reset_trackers()

            if t % self.params.outputs.print_frequency == 0:
                if self.params.outputs.network.calc_network_stats:
                    self.network_utils.write_network_stats(t=t)

                if self.params.outputs.network.calc_component_stats:
                    ao.print_components(
                        run_id,
                        t,
                        self.run_seed,
                        self.pop.pop_seed,
                        self.pop.connected_components(),
                        outdir,
                    )
                if self.params.outputs.network.edge_list:
                    path = os.path.join(outdir, "network", f"Edgelist_t{t}.txt")
                    self.network_utils.write_graph_edgelist(path)

        return stats

    def update_high_risk(self):
        """
        :Purpose:
            Update high risk agents or remove them from high risk pool
        """
        for agent in self.pop.high_risk_agents:
            if agent.high_risk_time > 0:
                agent.high_risk_time -= 1
                if (
                    agent.so == "HM"
                    and self.params.features.prep
                    and (
                        self.params.prep.target_model
                        in ("high_risk", "incarcerated_high_risk")
                    )
                ):
                    for part in agent.partners:
                        if not (part.hiv or part.vaccine):
                            self.initiate_prep(part, time)
            else:
                self.pop.high_risk_agents.remove_agent(agent)
                agent.high_risk = False

                if (
                    self.params.features.incar
                ):  # REVIEWED why does this check hm/hf and then subtract things - could this be more generic? Sarah to look into if this needs to be sex based
                    agent.mean_num_partners -= self.params.high_risk.partner_scale

    def initialize_random_trial(self):
        """
        :Purpose:
            Initialize random trial in population
        """
        assert params.model.network.enable, "Network must be enabled for random trial"

        print("Starting random trial")
        components = self.pop.connected_components()

        total_nodes = 0
        print(
            "Number of components",
            len([1 for comp in components if comp.number_of_nodes()]),
        )
        for comp in components:
            total_nodes += comp.number_of_nodes()

            if self.run_random.random() < 0.5:
                # Component selected as treatment pod!
                if not self.params.features.pca:
                    for ag in comp.nodes():
                        if not ag.hiv and not ag.prep:
                            ag.intervention_ever = True
                            if (
                                self.run_random.random() < self.params.prep.target
                                and not ag.vaccine
                            ):
                                self.initiate_prep(ag, time, force=True)
                elif self.params.prep.pca.choice == "eigenvector":
                    centrality = nx.algorithms.centrality.eigenvector_centrality(comp)
                    assert len(centrality) >= 1, "Empty centrality"
                    ordered_centrality = sorted(centrality, key=centrality.get)
                    intervention_agent = False
                    for ag in ordered_centrality:
                        if not ag.hiv:
                            ag.prep_awareness = True
                            ag.pca = True
                            ag.pca_suitable = True
                            intervention_agent = True
                            break
                    if not intervention_agent:
                        ag = ordered_centrality[0]
                        ag._pca = True
                elif self.params.prep.pca.choice == "bridge":
                    # list all edges that are bridges
                    all_bridges = list(nx.bridges(comp))
                    comp_agents = [
                        agent
                        for agents in all_bridges
                        for agent in agents
                        if not agent.hiv
                    ]  # all suitable agents in bridges
                    if comp_agents:
                        chosen_agent = self.run_random.choice(
                            comp_agents
                        )  # select change agent
                        chosen_agent.prep_awareness = True  # make aware
                        chosen_agent.pca = True
                        chosen_agent.pca_suitable = True
                    else:
                        chosen_agent = list(comp.nodes)[0]
                        chosen_agent.pca = True

                elif self.params.prep.pca.choice == "random":
                    suitable_agent_choices = [ag for ag in comp.nodes if not ag.hiv]
                    if (
                        suitable_agent_choices
                    ):  # if there are agents who meet eligibility criteria,
                        # select one randomly
                        chosen_agent = self.run_random.choice(suitable_agent_choices)
                        chosen_agent.pca = True
                        chosen_agent.pca_suitable = True
                        chosen_agent.prep_awareness = True  # make aware
                    else:  # if no suitable agents, mark a non-suitable agent
                        chosen_agent = self.run_random.choice(list(comp.nodes))
                        chosen_agent.pca = True

        print(("Total agents in trial: ", total_nodes))

    def update_all_agents(self, time: int, burn: bool = False):
        """
        :Purpose:
            Update PWID agents:
            For each agent:
                1 - determine agent type
                2 - get partners
                3 - agent interacts with partners
                5 - VCT (Voluntsry Counseling and Testing)
                6 - if PWID: needle exchange
                7 - if HIV: HAART, AIDS

        :Input:
            agent, time

        :Output:
            none
        """
        if time > 0 and self.params.features.static_n is False:
            self.pop.update_partner_assignments()

        for rel in self.pop.relationships:
            # If in burn, ignore interactions
            if not burn:
                self.agents_interact(time, rel)

            # If static network, ignore relationship progression
            if not self.params.features.static_n:
                if rel.progress():
                    self.pop.remove_relationship(rel)

        if (
            self.params.features.high_risk
        ):  # TO_REVIEW can this move into the general loop with a check for agent.high_risk? (agents should never become high risk if the feature isn't enabled)
            self.update_high_risk()

        for agent in self.pop.all_agents:
            if (
                self.params.features.pca
                and self.run_random.random() < self.params.prep.prep_awareness.prob
                and not burn
            ):
                agent.prep_awareness = True
                if self.run_random.random() < self.params.prep.pca.prep.prob:
                    self.initiate_prep(agent, time, force=True)

            if self.params.features.incar:
                self.incarcerate(agent, time)

            if agent.msmw and self.run_random.random() < self.params.msmw.hiv.prob:
                self.hiv_convert(agent)

            if agent.hiv:
                # If in burnin, ignore HIV
                if not burn:
                    self.diagnose_hiv(agent, time)
                    self.progress_to_aids(agent)

                    if self.params.features.haart:
                        self.update_haart(agent, time)
                        agent.hiv_time += 1
            else:
                if self.params.features.prep:
                    if time >= self.params.prep.start:
                        if agent.prep:
                            self.discontinue_prep(agent)
                        elif self.params.prep.target_model == "RandomTrial":
                            pass
                        elif agent.prep_eligible(self.params.prep.target_model):
                            self.initiate_prep(agent, time)
                    if self.params.features.vaccine and not agent.prep:
                        self.advance_vaccine(
                            agent, time, vaxType=self.params.vaccine.type, burn=burn
                        )

        if (
            self.params.features.prep
            and time == self.params.prep.start
            and self.params.prep.target_model == "RandomTrial"
        ):
            self.initialize_random_trial()

    def agents_interact(self, time: int, rel: Relationship):
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

        if self.params.features.pca:
            if rel.bond_type == "sexOnly" and rel.duration != rel.total_duration:
                pass
            else:
                self.pca_interaction(rel, time)

        # Agent 1 is HIV, partner is succept
        if rel.agent1.hiv and not rel.agent2.hiv:
            agent = rel.agent1
            partner = rel.agent2
        # If agent_2 is HIV agen1 is not, agent_2 is HIV, agent_1 is succept
        elif not rel.agent1.hiv and rel.agent2.hiv:
            agent = rel.agent2
            partner = rel.agent1
        else:  # neither agent is HIV or both are
            return False

        rel_sex_possible = sex_possible(agent.so, partner.so, self.params)
        partner_drug_type = partner.drug_use
        agent_drug_type = agent.drug_use

        if partner_drug_type == "Inj" and agent_drug_type == "Inj":
            # Injection is possible
            if rel_sex_possible:
                # Sex is possible
                rv = self.run_random.random()  # REVIEW after bond types established
                if rv < 0.25:  # Needle only (60%)
                    self.needle_transmission(agent, partner, time)
                else:  # Both sex and needle (20%)
                    self.needle_transmission(agent, partner, time)
                    self.sex_transmission(rel, time)
            else:
                # Sex not possible, needle only
                self.needle_transmission(agent, partner, time)

        elif partner_drug_type in ["NonInj", "None"] or agent_drug_type in [
            "NonInj",
            "None",
        ]:
            if rel_sex_possible:
                self.sex_transmission(rel, time)
            else:
                return False

        return True

    def pca_interaction(self, relationship: Relationship, time, force=False):
        """
        :Purpose:
            Simulate peer change agent interactions
            Knowledge if one agent is aware and one unaware, opinion if one agent swayint the other
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
            agent_opinion = agent.prep_opinion
            partner_opinion = partner.prep_opinion
            agent_influence = nx.closeness_centrality(self.pop.graph, agent)
            partner_influence = nx.closeness_centrality(self.pop.graph, partner)

            if agent_influence > partner_influence:
                partner.prep_opinion = np.mean(
                    [agent.prep_opinion, partner.prep_opinion]
                )
            elif agent_influence < partner_influence:
                agent.prep_opinion = np.mean([agent.prep_opinion, partner.prep_opinion])

            if self.run_random.random() < self.params.prep.pca.prep.prob:
                if (
                    agent_opinion
                    < self.params.prep.pca.prep_opinion.threshold
                    < agent.prep_opinion
                ):
                    self.initiate_prep(agent, time, force=True)
                elif (
                    partner_opinion
                    < self.params.prep.pca.prep_opinion.threshold
                    < partner.prep_opinion
                ):
                    self.initiate_prep(partner, time, force=True)

        def knowledge_dissemination(partner):
            partner.prep_awareness = True
            if (
                partner.prep_opinion > self.params.prep.pca.opinion.threshold
                and self.run_random.random() < self.params.prep.pca.prep.prob
            ):
                self.initiate_prep(partner, time, force=True)

        def transmission_probability():
            if (
                relationship.agent1.prep_awareness
                and relationship.agent2.prep_awareness
            ):
                p = self.params.prep.pca.opinion.transmission
            else:
                p = self.params.prep.pca.knowledge.transmission

            if num_acts == 1:
                p_total_transmission = p
            elif num_acts >= 1:
                p_total_transmission = 1.0 - binom.pmf(0, num_acts, p)
            else:
                p_total_transmission = 0

            return p_total_transmission

        acts_prob = self.run_random.random()
        acts_bin = 0
        current_p_value = 0.0

        while acts_prob > current_p_value:
            acts_bin += 1
            current_p_value += self.params.partnership.interaction[
                relationship.bond_type
            ][acts_bin].prob

        min = self.params.partnership.interaction[relationship.bond_type][acts_bin].min
        max = self.params.partnership.interaction[relationship.bond_type][acts_bin].max
        if min == max:
            num_acts = min
        else:
            num_acts = self.run_random.randrange(min, max)

        if num_acts < 1:
            return

        if (
            relationship.agent1.prep_awareness
            and not relationship.agent2.prep_awareness
        ):
            if self.run_random.random() < transmission_probability() or force:
                knowledge_dissemination(relationship.agent2)
        elif (
            not relationship.agent1.prep_awareness
            and relationship.agent2.prep_awareness
        ):
            if self.run_random.random() < transmission_probability() or force:
                knowledge_dissemination(relationship.agent1)
        elif (
            relationship.agent1.prep_awareness
            and relationship.agent2.prep_awareness
            or force
        ):
            if self.run_random.random() < transmission_probability() or force:
                influence(relationship.agent1, relationship.agent2)

    def needle_transmission(self, agent: Agent, partner: Agent, time: int):
        """
        :Purpose:
            Simulate random transmission of HIV between two PWID agents
            through needle.\n
            Agent must by HIV+ and partner not.

        :Input:
            agents : int
            partner : int
            time : int
        :Output: -
        """

        assert agent.hiv
        assert not partner.hiv
        assert agent.drug_use == "Inj"
        assert partner.drug_use == "Inj"

        agent_race = agent.race
        agent_sex_type = agent.so

        # REVIEWED why is the mean number of sex acts for a class multiplied by needle calibration? - change to num_needle_acts
        mean_num_acts = (
            self.params.demographics[agent_race][agent_sex_type].num_needle_acts
            * self.params.calibration.needle.act
        )
        share_acts = round(poisson.rvs(mean_num_acts, size=1)[0])

        if agent.sne:  # safe needle exchange - minimal sharing
            p_unsafe_needle_share = 0.02  # minimal needle sharing
        else:  # they do share a needle

            # If sharing, minimum of 1 share act
            if share_acts < 1:
                share_acts = 1

            p_unsafe_needle_share = (
                self.params.demographics[agent_race][agent_sex_type].needle_sharing
                * self.params.needle_exchange.prevalence
            )

        for n in range(share_acts):
            if self.run_random.random() > p_unsafe_needle_share:
                share_acts -= 1

        if share_acts >= 1.0:
            p = agent.get_transmission_probability("NEEDLE", self.params)

            p_total_transmission: float
            if share_acts == 1:
                p_total_transmission = p
            else:
                p_total_transmission = 1.0 - binom.pmf(0, share_acts, p)

            if self.run_random.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                self.hiv_convert(partner)

    def sex_transmission(self, rel: Relationship, time):
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
            agent.get_number_of_sex_acts(self.run_random, self.params)
            * self.params.calibration.sex.act
        )
        total_sex_acts = round(poisson.rvs(mean_sex_acts, size=1)[0])

        # Get condom usage
        if self.params.high_risk.condom_use_type == "Race":
            p_safe_sex = self.params.demographics[agent.race][agent.so].safe_sex
        else:
            p_safe_sex = prob.safe_sex(rel.total_sex_acts)

        # Reduction of risk acts between partners for condom usage
        unsafe_sex_acts = total_sex_acts
        for n in range(unsafe_sex_acts):
            if self.run_random.random() < p_safe_sex:
                unsafe_sex_acts -= 1

        if unsafe_sex_acts >= 1:
            # agent is HIV+
            rel.total_sex_acts += unsafe_sex_acts
            p_per_act = agent.get_transmission_probability("SEX", self.params)

            # Reduction of transmissibility for acts between partners for PrEP adherence
            if agent.prep or partner.prep:
                if "Oral" in self.params.prep.type:  # params.prep.type == "Oral":
                    if agent.prep_adherence == 1 or partner.prep_adherence == 1:
                        p_per_act *= 1.0 - self.params.prep.efficacy.adherent  # 0.04
                    else:
                        p_per_act *= (
                            1.0 - self.params.prep.efficacy.non_adherant
                        )  # 0.24

                elif "Inj" in self.params.prep.type:
                    p_per_act_reduction = (
                        -1.0 * np.exp(-5.528636721 * partner.prep_load) + 1
                    )
                    if agent.prep_adherence == 1 or partner.prep_adherence == 1:
                        p_per_act *= 1.0 - p_per_act_reduction  # 0.04

            if partner.vaccine:
                p_per_act_perc: float
                if self.params.vaccine.type == "HVTN702":
                    p_per_act_perc *= np.exp(
                        -2.88 + 0.76 * (np.log((partner.vaccine_time + 0.001) * 30))
                    )
                elif self.params.vaccine.type == "RV144":
                    p_per_act_perc *= np.exp(
                        -2.40 + 0.76 * (np.log(partner.vaccine_time))
                    )

                p_per_act *= 1 - p_per_act_perc

            p_total_transmission: float
            if unsafe_sex_acts == 1:
                p_total_transmission = p_per_act
            else:
                p_total_transmission = 1.0 - binom.pmf(0, unsafe_sex_acts, p_per_act)

            if self.run_random.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                self.hiv_convert(partner)

    def hiv_convert(self, agent: Agent):  # TODO rename
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

    def enroll_needle_exchange(self):
        """
        :Purpose:
            Enroll PWID agents in needle exchange
        """
        print(("\n\n!!!!Engaginge treatment process"))
        self.needle_exchange = True
        for agent in self.pop.all_agents:
            if (
                self.run_random.random() < self.params.needle_exchange.coverage
                and agent.drug_use == "Inj"
            ):
                agent.sne = True
                agent.intervention_ever = True

    def become_high_risk(self, agent: Agent, duration: int = None):

        if not self.params.features.high_risk:
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
            agent.high_risk_time = self.params.high_risk.sex_based[agent.so].duration

    def incarcerate(self, agent: Agent, time: int):
        """
        :Purpose:
            To incarcerate an agent or update their incarceration variables

        :Input:
            agent : int
            time : int

        """
        if not self.params.features.incar:
            return None

        hiv_bool = agent.hiv

        if agent.incar:
            agent.incar_time -= 1

            if agent.incar_time == 0:  # FREE AGENT
                self.pop.incarcerated_agents.remove_agent(agent)
                self.new_incar_release.add_agent(agent)
                agent.incar = False
                agent.incar_ever = True
                if (
                    not agent.high_risk and self.params.features.high_risk
                ):  # If behavioral treatment on and agent HIV, ignore HR period.
                    self.become_high_risk(agent)

                    agent.mean_num_partners = (
                        agent.mean_num_partners + self.params.high_risk.partner_scale
                    )

                if hiv_bool:
                    if agent.haart:
                        if (
                            self.run_random.random()
                            <= self.params.incar.haart.discontinue
                        ):  # 12% remain surpressed
                            agent.haart = False
                            agent.haart_adherence = 0
                            self.pop.intervention_haart_agents.remove_agent(agent)

                        # END FORCE

        elif self.run_random.random() < (
            self.params.demographics[agent.race][agent.so].incar.prob
            * (1 + (hiv_bool * 4))
            * self.params.calibration.incarceration
        ):
            jail_duration = self.params.demographics[agent.race][
                agent.so
            ].incar.duration.prob

            bin = current_p_value = 1
            p = self.run_random.random()
            while p >= current_p_value:
                current_p_value += jail_duration[bin].prob
                bin += 1

            timestay = self.run_random.randint(
                jail_duration[bin].min, jail_duration[bin].max
            )

            if hiv_bool:
                if not agent.hiv_dx:
                    if self.run_random.random() < self.params.incar.hiv.dx:
                        agent.hiv_dx = True
                else:  # Then tested and HIV, check to enroll in ART
                    if self.run_random.random() < self.params.incar.haart.prob:
                        tmp_rnd = self.run_random.random()
                        haart_adh = self.params.incar.haart.adherence
                        if tmp_rnd < haart_adh:
                            adherence = 5
                        else:
                            adherence = self.run_random.randint(1, 4)

                        # Add agent to HAART class set, update agent params
                        agent.haart = True
                        agent.intervention_ever = True
                        agent.haart_adherence = adherence
                        agent.haart_time = time
                        self.pop.intervention_haart_agents.add_agent(agent)

            agent.incar = True
            agent.incar_time = timestay
            self.pop.incarcerated_agents.add_agent(agent)

            # PUT PARTNERS IN HIGH RISK
            for partner in agent.partners:
                if not partner.high_risk and self.params.features.high_risk:
                    if self.run_random.random() < self.params.high_risk.proportion:
                        self.become_high_risk(partner)

                if self.params.features.prep and (
                    self.params.prep.target_model in ("Incar", "IncarHR")
                ):
                    # Atempt to put partner on prep if less than probability
                    if not partner.hiv and not agent.vaccine:
                        self.initiate_prep(partner, time)

    # REVIEW - change verbage to diagnosed
    def diagnose_hiv(self, agent: Agent, time: int):
        """
        :Purpose:
            Test the agent for HIV. If detected, add to identified list.

        :Input:
            agent : agent_Class
            time : int

        :Output:
            none
        """
        sex_type = agent.so
        race_type = agent.race
        tested = agent.hiv_dx

        def diagnose(agent):
            agent.hiv_dx = True
            self.new_dx.add_agent(agent)
            self.pop.intervention_dx_agents.add_agent(agent)
            if (
                self.params.features.partner_tracing
            ):  # TODO fix this logic; should get partnerTraced and then lose it after
                # For each partner, determine if found by partner testing
                for ptnr in agent.partners:
                    if ptnr.hiv and not ptnr.hiv_dx:
                        ptnr.partner_traced = True
                        ptnr.trace_time = time + 1

        if not tested:
            test_prob = self.params.demographics[race_type][sex_type].hiv.dx.prob

            # Rescale based on calibration param
            test_prob *= self.params.calibration.test_frequency

            # If roll less than test probablity
            if self.run_random.random() < test_prob:
                # Become tested, add to tested agent set
                diagnose(agent)
                # If treatment co-enrollment enabled and coverage greater than 0

            elif (
                agent.partner_traced
                and self.run_random.random() < 0.87
                and agent.trace_time == time
            ):
                diagnose(agent)

        agent.partner_traced = False

    def update_haart(self, agent: Agent, time: int):
        """
        :Purpose:
            Account for HIV treatment through highly active antiretroviral therapy (HAART).
            HAART was implemented in 1996, hence, there is treatment only after 1996.
            HIV treatment assumes that the agent knows their HIV+ status.

        :Input:
            agent : Agent
            time : int

        :Output:
            none
        """
        if not self.params.features.haart:
            return None

        # Check valid input
        assert agent.hiv

        agent_haart = agent.haart
        agent_race = agent.race
        agent_so = agent.so

        # Determine probability of HIV treatment
        if time >= 0 and agent.hiv_dx:
            # Go on HAART
            if not agent_haart and agent.haart_time == 0:
                if self.run_random.random() < (
                    self.params.demographics[agent_race][agent_so].haart.prob
                    * self.params.calibration.haart_coverage
                ):

                    haart_adh = self.params.demographics[agent_race][
                        agent_so
                    ].haart.adherence
                    if self.run_random.random() < haart_adh:
                        adherence = 5
                    else:
                        adherence = self.run_random.randint(1, 4)

                    # Add agent to HAART class set, update agent params
                    agent.haart = True
                    agent.intervention_ever = True
                    agent.haart_adherence = adherence
                    agent.haart_time = time
                    self.pop.intervention_haart_agents.add_agent(agent)

            # Go off HAART
            elif (
                agent_haart
                and self.run_random.random()
                < self.params.demographics[agent_race][agent_so].haart.discontinue
            ):
                agent.haart = False
                agent.haart_adherence = 0
                agent.haart_time = 0
                self.pop.intervention_haart_agents.remove_agent(agent)

    def discontinue_prep(self, agent: Agent, force: bool = False):
        # Agent must be on PrEP to discontinue PrEP
        assert agent.prep

        # If force flag set, auto kick off prep.
        if force:
            self.pop.intervention_prep_agents.remove_agent(agent)
            self.prep_agents[agent.race][agent.so] -= 1
            agent.prep = False
            agent.prep_reason = []

        # else if agent is on PrEP, see if they should discontinue
        else:
            if (
                self.run_random.random()
                < self.params.demographics[agent.race][agent.so].prep.discontinue
            ):
                self.pop.intervention_prep_agents.remove_agent(agent)
                self.prep_agents[agent.race][agent.so] -= 1

                if "Oral" in self.params.prep.type:
                    agent.prep = False
                    agent.prep_type = ""
                    agent.prep_reason = []
            else:  # if not discontinue, see if its time for a new shot. # REVIEWED what is this logic doing? This decrements, then update_prep_load increments - sarah to review with max
                if agent.prep_last_dose > 2:
                    agent.prep_last_dose = -1

        if self.params.prep.type == "Inj":
            agent.update_prep_load(self.params)

    def advance_vaccine(self, agent: Agent, time: int, vaxType: str, burn: bool):
        """
        :Purpose:
            Progress vaccine. Agents may receive injection or progress in time since injection.

        :Input:
            agent: Agent
            time: int

        :Output:
            none
        """
        if not self.params.features.vaccine:
            return None

        if agent.vaccine and not burn:
            agent.vaccine_time += 1
            if (
                self.params.vaccine.booster
                and agent.vaccine_time
                == self.params.demographics[agent.race][
                    agent.so
                ].vaccine.booser.interval
                and self.run_random.random()
                < self.params.demographics[agent.race][agent.so].vaccine.booster.prob
            ):
                agent.vaccinate(vaxType)

        elif time == self.params.vaccine.start:
            if self.params.vaccine.init == burn:  # both true or both false
                if (
                    self.run_random.random()
                    < self.params.demographics[agent.race][agent.so].vaccine.prob
                ):
                    agent.vaccinate(vaxType)

    def initiate_prep(self, agent: Agent, time: int, force: bool = False):
        """
        :Purpose:
            Place agents onto PrEP treatment.
            PrEP treatment assumes that the agent knows their HIV+ status is negative.

        :Input:
            agent : Agent
            time : int
            force : default is `False`

        :Output:
            none
        """

        def enroll_prep(self, agent: Agent):
            agent.prep = True
            agent.intervention_ever = True
            self.pop.intervention_prep_agents.add_agent(agent)
            self.new_prep.add_agent(agent)

            self.prep_agents[agent.race][agent.so] += 1

            if (
                self.run_random.random()
                < self.params.demographics[agent.race][agent.so].prep.adherence
            ):
                agent.prep_adherence = 1
            else:
                agent.prep_adherence = 0

            # set PrEP load and dosestep for PCK
            if "Inj" in self.params.prep.type and "Oral" in self.params.prep.type:
                agent.prep_load = self.params.prep.peak_load
                agent.prep_last_dose = 0

                if self.run_random.random() < self.params.prep.lai.prob:
                    agent.prep_type = "Inj"
                else:
                    agent.prep_type = "Oral"

        # agent must exist
        assert agent is not None

        # Prep only valid for agents not on prep and are HIV negative
        if agent.prep or agent.hiv:
            return

        # Determine probability of HIV treatment
        if force:
            enroll_prep(self, agent)
        else:
            if self.params.prep.target_model == "Racial":
                num_prep_agents = sum(self.prep_agents[agent.race].values())
            else:
                num_prep_agents = self.pop.intervention_prep_agents.num_members()

            if self.params.prep.target_model in ("Incar", "IncarHR"):
                if self.run_random.random() < self.params.prep.target:
                    enroll_prep(self, agent)
                return None
            elif self.params.prep.target_model == "Racial":
                all_hiv_agents = set(self.pop.all_agents.subset["HIV"].members)
                all_race = set(
                    self.pop.all_agents.subset["Race"].subset[agent.race].members
                )
                hiv_agents = len(all_hiv_agents & all_race)
                target_prep = (
                    int(
                        self.pop.all_agents.subset["Race"]
                        .subset[agent.race]
                        .num_members()
                    )
                    - hiv_agents
                ) * self.params.demographics[agent.race][agent.so].prep.coverage

            else:
                target_prep = int(
                    (
                        self.pop.all_agents.num_members()
                        - self.pop.all_agents.subset["HIV"].num_members()
                    )
                    * self.params.prep.target
                )

            if self.params.prep.target_model in ("Incar", "IncarHR"):
                if self.run_random.random() < self.params.prep.target:
                    enroll_prep(self, agent)
            elif (
                num_prep_agents < target_prep
                and time >= self.params.prep.start
                and agent.prep_eligible(self.params.prep.target_model)
            ):
                enroll_prep(self, agent)

    def progress_to_aids(self, agent: Agent):
        """
        :Purpose:
            Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        assert agent.hiv

        # REVIEWED Why do we check for not HAART, but then get HAART adherance? - Sarah to ask Max
        if not agent.haart:
            p = prob.adherence_prob(agent.haart_adherence)

            if self.run_random.random() < p * self.params.hiv.aids.prob:
                agent.aids = True
                self.pop.hiv_aids_agents.add_agent(agent)

    def die_and_replace(self):

        """
        :Purpose:
            Let agents die and replace the dead agent with a new agent randomly.
        """
        # die stage
        for agent in self.pop.all_agents.members:

            # agent incarcerated, don't evaluate for death
            if agent.incar:
                continue

            # death rate per 1 person-month
            p = (
                prob.get_death_rate(
                    agent.hiv,
                    agent.aids,
                    agent.race,
                    agent.haart_adherence,
                    self.params,
                )
                * self.params.calibration.mortality
            )

            if self.run_random.random() < p:
                self.deaths.append(agent)

                # End all existing relationships
                for rel in agent.relationships:
                    rel.progress(force=True)
                    self.pop.remove_relationship()

        # replace stage
        for agent in self.deaths:
            # Remove agent from agent class and sub-sets
            self.pop.remove_agent(agent)

            new_agent = self.pop.create_agent(agent.race, agent.so)
            self.pop.add_agent(new_agent)
