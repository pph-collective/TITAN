#!/usr/bin/env python3
# encoding: utf-8

# Imports
import random
from random import Random
from typing import Dict, Any, List, Sequence, Optional
import os
import uuid

import numpy as np  # type: ignore
from scipy.stats import binom  # type: ignore
from scipy.stats import poisson  # type: ignore
import networkx as nx  # type: ignore
from dotmap import DotMap  # type: ignore


from .agent import AgentSet, Agent, Relationship
from .network_graph_tools import NetworkClass
from . import analysis_output as ao
from . import probabilities as prob
from .ABM_partnering import sex_possible


class HIVModel(NetworkClass):
    """
    :Purpose:
        This is the core class used to simulate
        the spread of HIV and drug use in one MSA
        (Metropolitan Statistical Area).

    :Input:
        params: DotMap - the parameter object for this model
    """

    def __repr__(self):
        returnStr = "\n"
        returnStr += "Seed: %d\n" % (self.runseed)
        returnStr += "Npop: %d\n" % (self.params.model.num_pop)
        returnStr += "Time: %d\n" % (self.params.model.time_range)

        return returnStr

    def __init__(self, params: DotMap):

        self.params = params

        def get_check_rand_int(seed):
            """
            Check the value passed of a seed, make sure it's an int, if 0, get a random seed
            """
            if type(seed) is not int:
                raise ValueError("Random seed must be integer")
            elif seed == 0:
                return random.randint(1, 1000000)
            else:
                return seed

        self.runseed = get_check_rand_int(params.model.seed.run)
        self.popseed = get_check_rand_int(params.model.seed.ppl)
        self.netseed = get_check_rand_int(params.model.seed.net)

        print("=== Begin Initialization Protocol ===\n")

        NetworkClass.__init__(self, params, popSeed=self.popseed, netSeed=self.netseed)

        print("\n\tCreating lists")
        # Other lists / dictionaries
        self.NewInfections = AgentSet("NewInfections")
        self.NewDiagnosis = AgentSet("NewDiagnosis")
        self.NewIncarRelease = AgentSet("NewIncarRelease")
        self.NewHRrolls = AgentSet("NewHRrolls")
        self.newPrEPagents = AgentSet("NewPrEPagents")

        self.totalDiagnosis = 0
        self.needle_exchange = False

        self.prep_agents: Dict[str, Dict[str, int]] = {}
        for race in params.classes.races:
            self.prep_agents[race] = {}
            for st in params.classes.sex_types:
                self.prep_agents[race][st] = 0

        # Set seed format. 0: pure random, -1: Stepwise from 1 to nRuns, else: fixed value
        print(("\tRun seed was set to:", self.runseed))
        self.runRandom = Random(self.runseed)
        random.seed(self.runseed)
        np.random.seed(self.runseed)
        print(("\tFIRST RANDOM CALL %d" % random.randint(0, 100)))

        print("\tResetting death count")
        self.deathSet: List[Agent] = []  # Number of death

        print("\tCreating network graph")
        self.create_graph_from_agents(
            self.All_agentSet
        )  # REVIEWED redundant with NetworkClass init? - review with max, logic feels scattered as NetworkClass also intializes a graph

        print("\n === Initialization Protocol Finished ===")

    def run(self, outdir):
        """
        Core of the model:
            1. Prints networkReport for first agents.
            2. Makes agents become HIV (used for current key_time tracking for acute)
            3. Loops over all time steps
                a. _update AllAgents()
                b. reset death count
                c. _ self._die_and_replace()
                d. self._update_population()
                e. self._reset_partner_count()
        """

        def print_stats(stat: Dict[str, Dict[str, int]], run_id: uuid.UUID):
            for report in self.params.outputs.reports:
                printer = getattr(ao, report)
                printer(
                    run_id,
                    t,
                    self.runseed,
                    self.popseed,
                    self.netseed,
                    stat,
                    self.params,
                    outdir,
                )

        def get_components():
            return sorted(
                nx.connected_components(self.get_Graph()), key=len, reverse=True
            )

        def reset_trackers():
            self.NewInfections.clear_set()
            self.NewDiagnosis.clear_set()
            self.NewHRrolls.clear_set()
            self.NewIncarRelease.clear_set()
            self.newPrEPagents.clear_set()
            self.deathSet = []

        def burnSimulation(burnDuration: int):
            print(
                ("\n === Burn Initiated for {} timesteps ===".format(burnDuration + 1))
            )
            for t in range(0, burnDuration + 1):
                self._update_AllAgents(t, burn=True)

                if self.params.features.die_and_replace:
                    self._die_and_replace()
            print(("\tBurn Cuml Inc:\t{}".format(self.NewInfections.num_members())))
            reset_trackers()
            print(" === Simulation Burn Complete ===")

        def makeAgentZero(numPartners: int):
            firstHIV = self.runRandom.choice(self.DU_Inj_agentSet.members)
            for i in range(numPartners):
                self.update_agent_partners(self.G, firstHIV, self.params)
            self._become_HIV(firstHIV)

        run_id = uuid.uuid4()

        burnSimulation(self.params.model.burn_duration)

        print("\n === Begin Simulation Run ===")
        if self.params.outputs.draw_figures:
            nNodes = self.G.number_of_nodes()
            self.visualize_network(
                node_size=5000.0 / nNodes,
                curtime=0,
                iterations=10,
                label="Seed" + str(self.runseed),
            )

        if self.params.outputs.calc_component_stats:
            ao.print_components(
                run_id,
                0,
                self.runseed,
                self.popseed,
                self.netseed,
                get_components(),
                outdir,
            )

        print("\t===! Start Main Loop !===")

        # dictionary to hold results over time
        stats = {}

        # If we are using an agent zero method, create agent zero.
        if self.params.features.agent_zero:
            makeAgentZero(4)

        if self.params.outputs.edge_list:
            path = os.path.join(outdir, "network", f"Edgelist_t0.txt")
            self.write_G_edgelist(path)

        for t in range(1, self.params.model.time_range + 1):
            print(f"\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t.: TIME {t}")
            if (
                self.params.outputs.draw_figures
                and t % self.params.outputs.print_frequency == 0
            ):
                self.visualize_network(
                    node_size=5000.0 / nNodes,
                    curtime=t,
                    iterations=10,
                    label="Seed" + str(self.runseed),
                )
            # todo: GET THIS TO THE NEW HIV COUNT

            print(
                "\tSTARTING HIV count:{}\tTotal Incarcerated:{}\tHR+:{}\tPrEP:{}".format(
                    self.HIV_agentSet.num_members(),
                    self.incarcerated_agentSet.num_members(),
                    self.highrisk_agentsSet.num_members(),
                    self.Trt_PrEP_agentSet.num_members(),
                )
            )

            self._update_AllAgents(t)

            if self.params.features.die_and_replace:
                self._die_and_replace()

            stats[t] = ao.get_stats(
                self.All_agentSet,
                self.HIV_agentSet,
                self.incarcerated_agentSet,
                self.Trt_PrEP_agentSet,
                self.newPrEPagents,
                self.NewInfections,
                self.NewDiagnosis,
                self.Relationships,
                self.NewHRrolls,
                self.NewIncarRelease,
                self.deathSet,
                self.params,
            )
            print_stats(stats[t], run_id)

            print(("Number of relationships: %d" % len(self.Relationships)))
            self.All_agentSet.print_subsets()

            self.totalDiagnosis += len(self.NewDiagnosis.members)
            if (
                self.totalDiagnosis > self.params.needle_exchange.init_at_pop
                and not self.needle_exchange
            ):
                self.enroll_needle_exchange()

            # RESET counters for the next time step
            reset_trackers()

            if t % self.params.outputs.print_frequency == 0:
                if self.params.outputs.calc_network_stats:
                    self.write_network_stats(t=t)
                if self.params.outputs.calc_component_stats:
                    ao.print_components(
                        run_id,
                        t,
                        self.runseed,
                        self.popseed,
                        self.netseed,
                        get_components(),
                        outdir,
                    )
                if self.params.outputs.edge_list:
                    path = os.path.join(outdir, "network", f"Edgelist_t{t}.txt")
                    self.write_G_edgelist(path)

        return stats

    def _update_AllAgents(self, time: int, burn: bool = False):
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
            self.update_partner_assignments(self.G, self.params)

        for rel in self.Relationships:
            # If in burn, ignore interactions
            if not burn:
                self._agents_interact(rel)

            # If static network, ignore relationship progression
            if not self.params.features.static_n:
                if rel.progress():
                    g = self.G
                    if g.has_edge(rel.agent1, rel.agent2):
                        g.remove_edge(rel.agent1, rel.agent2)

                    self.Relationships.remove(rel)
                    del rel

        if self.params.features.high_risk:
            for agent in self.highrisk_agentsSet:
                if agent.high_risk_time > 0:
                    agent.high_risk_time -= 1
                    if (
                        agent.so == "HM"
                        and self.params.features.prep
                        and (self.params.prep.target_model in ("HR", "IncarHR"))
                    ):
                        for part in agent.partners:
                            if not (part.hiv or part.vaccine):
                                self._initiate_PrEP(part, time)
                else:
                    self.highrisk_agentsSet.remove_agent(agent)
                    agent.high_risk = False

                    if (
                        self.params.features.incar
                    ):  # REVIEWED why does this check hm/hf and then subtract things - could this be more generic? Sarah to look into if this needs to be sex based
                        agent.neam_num_partners -= self.params.high_risk.partner_scale

        for agent in self.All_agentSet:
            if self.params.features.incar:
                self._incarcerate(agent, time)

            if agent.msmw and self.runRandom.random() < self.params.msmw.hiv.prob:
                self._become_HIV(agent)

            if agent.hiv:
                # If in burnin, ignore HIV
                if not burn:
                    self._HIVtest(agent, time)
                    self._progress_to_AIDS(agent)

                    if self.params.features.haart:
                        self._update_HAART(agent, time)
                        agent.hiv_time += 1
            else:
                if self.params.features.prep:
                    if time >= self.params.prep.start:
                        if agent.prep:
                            self._discont_PrEP(agent)
                        elif self.params.prep.target_model == "RandomTrial":
                            pass
                        elif agent.PrEP_eligible(self.params.prep.target_model):
                            self._initiate_PrEP(agent, time)
                    if self.params.features.vaccine and not agent.prep and not burn:
                        self.advance_vaccine(
                            agent, time, vaxType=self.params.vaccine.type
                        )

        if self.params.features.prep and time >= self.params.prep.start:
            if (
                self.params.prep.target_model == "RandomTrial"
                and time == self.params.prep.start
            ):
                print("Starting random trial")
                components = sorted(self.connected_components(), key=len, reverse=True)
                totNods = 0
                for comp in components:
                    totNods += comp.number_of_nodes()
                    if self.runRandom.random() < 0.5:
                        # Component selected as treatment pod!
                        for ag in comp.nodes():
                            if (not ag.hiv) and (not ag.prep):
                                ag.intervention_ever = True
                                if (
                                    self.runRandom.random() < self.params.prep.target
                                    and not agent.vaccine
                                ):
                                    self._initiate_PrEP(ag, time, force=True)
                print(("Total agents in trial: ", totNods))

    def _agents_interact(self, rel: Relationship) -> bool:
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

        if rel.agent1.hiv and not rel.agent2.hiv:  # Agent 1 is HIV, partner is succept
            agent = rel.agent1
            partner = rel.agent2
        elif (
            not rel.agent1.hiv and rel.agent2.hiv
        ):  # If agent_2 is HIV agen1 is not, agent_2 is HIV, agent_1 is succept
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
                rv = self.runRandom.random()
                if rv < 0.25:  # Needle only (60%)
                    self._needle_transmission(agent, partner)
                else:  # Both sex and needle (20%)
                    self._needle_transmission(agent, partner)
                    self._sex_transmission(rel)
            else:
                # Sex not possible, needle only
                self._needle_transmission(agent, partner)

        elif partner_drug_type in ["NonInj", "None"] or agent_drug_type in [
            "NonInj",
            "None",
        ]:
            if rel_sex_possible:
                self._sex_transmission(rel)
            else:
                return False

        return True  # if didn't short circuit, agents interacted

    def _needle_transmission(self, agent: Agent, partner: Agent):
        """
        :Purpose:
            Simulate random transmission of HIV between two PWID agents
            through needle.\n
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

        # REVIEWED why is the mean number of sex acts for a class multiplied by needle calibration? - change to num_needle_acts
        MEAN_N_ACTS = (
            self.params.demographics[agent_race][agent_sex_type].num_needle_acts
            * self.params.calibration.needle.act
        )
        share_acts = round(poisson.rvs(MEAN_N_ACTS, size=1)[0])

        if agent.sne:  # safe needle exchange - minimal sharing
            p_UnsafeNeedleShare = 0.02  # minimal needle sharing
        else:  # they do share a needle

            # If sharing, minimum of 1 share act
            if share_acts < 1:
                share_acts = 1

            p_UnsafeNeedleShare = (
                self.params.demographics[agent_race][agent_sex_type].needle_sharing
                * self.params.needle_exchange.prevalence
            )

        for n in range(share_acts):
            if self.runRandom.random() > p_UnsafeNeedleShare:
                share_acts -= 1

        if share_acts >= 1.0:
            p = agent.get_transmission_probability("NEEDLE", self.params)

            p_total_transmission: float
            if share_acts == 1:
                p_total_transmission = p
            else:
                p_total_transmission = 1.0 - binom.pmf(0, share_acts, p)

            if self.runRandom.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                self._become_HIV(partner)

    def _sex_transmission(self, rel: Relationship):
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
        MSexActs = (
            agent.get_number_of_sex_acts(self.runRandom, self.params)
            * self.params.calibration.sex.act
        )
        T_sex_acts = round(poisson.rvs(MSexActs, size=1)[0])

        # Get condom usage
        if self.params.high_risk.condom_use_type == "Race":
            p_SafeSex = self.params.demographics[agent.race][agent.so].safe_sex
        else:
            p_SafeSex = prob.safe_sex(rel.total_sex_acts)

        # Reduction of risk acts between partners for condom usage
        U_sex_acts = T_sex_acts
        for n in range(U_sex_acts):
            if self.runRandom.random() < p_SafeSex:
                U_sex_acts -= 1

        if U_sex_acts >= 1:
            # agent is HIV+
            rel.total_sex_acts += U_sex_acts
            ppAct = agent.get_transmission_probability("SEX", self.params)

            # Reduction of transmissibility for acts between partners for PrEP adherence
            if agent.prep or partner.prep:
                if "Oral" in self.params.prep.type:  # params.prep.type == "Oral":
                    if agent.prep_adherence == 1 or partner.prep_adherence == 1:
                        ppAct = ppAct * (
                            1.0 - self.params.prep.efficacy.adherent
                        )  # 0.04
                    else:
                        ppAct = ppAct * (
                            1.0 - self.params.prep.efficacy.non_adherant
                        )  # 0.24

                elif "Inj" in self.params.prep.type:
                    ppActReduction = -1.0 * np.exp(-5.528636721 * partner.prep_load) + 1
                    if agent.prep_adherence == 1 or partner.prep_adherence == 1:
                        ppAct = ppAct * (1.0 - ppActReduction)  # 0.04

            if partner.vaccine:
                if self.params.vaccine.type == "HVTN702":
                    ppAct *= np.exp(
                        -2.88 + 0.76 * (np.log((partner.vaccine_time + 0.001) * 30))
                    )

                elif self.params.vaccine.type == "RV144":
                    ppAct *= np.exp(-2.40 + 0.76 * (np.log(partner.vaccine_time)))

            p_total_transmission: float
            if U_sex_acts == 1:
                p_total_transmission = ppAct
            else:
                p_total_transmission = 1.0 - binom.pmf(0, U_sex_acts, ppAct)

            if self.runRandom.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                self._become_HIV(partner)

    def _become_HIV(self, agent: Agent):
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
            self.NewInfections.add_agent(agent)
            self.HIV_agentSet.add_agent(agent)

        if agent.prep:
            self._discont_PrEP(agent, force=True)

    def enroll_needle_exchange(self):
        """
        :Purpose:
            Enroll PWID agents in needle exchange
        """
        print(("\n\n!!!!Engaginge treatment process"))
        self.needle_exchange = True
        for agent in self.All_agentSet:
            if (
                self.runRandom.random() < self.params.needle_exchange.coverage
                and agent.drug_use == "Inj"
            ):
                agent.sne = True
                agent.intervention_ever = True

    # REVIEWED this isn't used anywhere, but should be! _incarcerate makes things high risk and should reference this
    def _becomeHighRisk(self, agent: Agent, duration: int = None):

        if agent not in self.highrisk_agentsSet.members:
            self.highrisk_agentsSet.add_agent(agent)

        if not agent.high_risk_ever:
            self.NewHRrolls.add_agent(agent)

        agent.high_risk = True
        agent.high_risk_ever = True

        if duration is not None:
            agent.high_risk_time = duration
        else:
            agent.high_risk_time = self.params.high_risk.sex_based[agent.so].duration

    def _incarcerate(self, agent: Agent, time: int):
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
                self.incarcerated_agentSet.remove_agent(agent)
                self.NewIncarRelease.add_agent(agent)
                agent.incar = False
                agent.incar_ever = True
                if (
                    not agent.high_risk and self.params.features.high_risk
                ):  # If behavioral treatment on and agent HIV, ignore HR period.
                    if (
                        self.params.incar.treatment.high_risk
                        and hiv_bool
                        and (time >= self.params.incar.treatment.start)
                    ):
                        pass
                    else:  # Else, become high risk
                        self.highrisk_agentsSet.add_agent(agent)
                        if not agent.high_risk_ever:
                            self.NewHRrolls.add_agent(agent)

                        agent.neam_num_partners = (
                            agent.neam_num_partners
                            + self.params.high_risk.partner_scale
                        )
                        agent.high_risk = True
                        agent.high_risk_ever = True
                        agent.high_risk_time = self.params.HR_M_dur

                if hiv_bool:
                    if agent.haart:
                        if (
                            self.runRandom.random()
                            > self.params.incar.haart.discontinue
                        ):  # 12% remain surpressed
                            pass

                        else:
                            agent.haart = False
                            agent.haart_adherence = 0
                            self.Trt_ART_agentSet.remove_agent(agent)

                        # END FORCE

        elif self.runRandom.random() < (
            self.params.demographics[agent.race][agent.so].incar.prob
            * (1 + (hiv_bool * 4))
            * self.params.calibration.incarceration
        ):
            # REVIEWED what about other sex types? -needs to be generalized - Sarah meeting with someone
            if agent.so == "HF":
                jailDuration = prob.HF_jail_duration
            elif agent.so == "HM":
                jailDuration = prob.HM_jail_duration

            durationBin = current_p_value = 1
            p = self.runRandom.random()
            while p >= current_p_value:
                current_p_value += jailDuration[durationBin]["p_value"]
                durationBin += 1

            timestay = self.runRandom.randint(
                jailDuration[durationBin]["min"], jailDuration[durationBin]["max"]
            )

            if hiv_bool:
                if not agent.hiv_dx:
                    if self.runRandom.random() < self.params.incar.hiv.dx:
                        agent.hiv_dx = True
                else:  # Then tested and HIV, check to enroll in ART
                    if self.runRandom.random() < self.params.incar.haart.prob:
                        tmp_rnd = self.runRandom.random()
                        haart_adh = self.params.incar.haart.adherence
                        if tmp_rnd < haart_adh:
                            adherence = 5
                        else:
                            adherence = self.runRandom.randint(1, 4)

                        # Add agent to HAART class set, update agent params
                        agent.haart = True
                        agent.intervention_ever = True
                        agent.haart_adherence = adherence
                        agent.haart_time = time
                        self.Trt_ART_agentSet.add_agent(agent)

            agent.incar = True
            agent.incar_time = timestay
            self.incarcerated_agentSet.add_agent(agent)

            # PUT PARTNERS IN HIGH RISK
            for partner in agent.partners:
                if not partner.high_risk:
                    if self.runRandom.random() < self.params.high_risk.proportion:
                        if not partner.high_risk:
                            self.highrisk_agentsSet.add_agent(partner)
                            if not partner.high_risk_ever:
                                self.NewHRrolls.add_agent(partner)
                            partner.neam_num_partners += (
                                self.params.high_risk.partner_scale
                            )  # 32.5 #2 + 3.25 from incar HR
                            partner.high_risk = True
                            partner.high_risk_ever = True
                            partner.high_risk_time = self.params.high_risk.sex_based[
                                partner.so
                            ].duration
                if self.params.features.prep and (
                    self.params.prep.target_model in ("Incar", "IncarHR")
                ):
                    # Atempt to put partner on prep if less than probability
                    if not partner.hiv and not agent.vaccine:
                        self._initiate_PrEP(partner, time)

    # REVIEW - change verbage to diagnosed
    def _HIVtest(self, agent: Agent, time: int):
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
            self.NewDiagnosis.add_agent(agent)
            self.Trt_Tstd_agentSet.add_agent(agent)
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
            if self.runRandom.random() < test_prob:
                # Become tested, add to tested agent set
                diagnose(agent)
                # If treatment co-enrollment enabled and coverage greater than 0

            elif (
                agent.partner_traced
                and self.runRandom.random() < 0.87
                and agent.trace_time == time
            ):
                diagnose(agent)

        agent.partner_traced = False

    def _update_HAART(self, agent: Agent, time: int):
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

        # Check valid input
        assert agent.hiv

        agent_haart = agent.haart
        agent_race = agent.race
        agent_so = agent.so

        # Determine probability of HIV treatment
        if time >= 0 and agent.hiv_dx:
            # Go on HAART
            if not agent_haart and agent.haart_time == 0:
                if self.runRandom.random() < (
                    self.params.demographics[agent_race][agent_so].haart.prob
                    * self.params.calibration.haart_coverage
                ):

                    haart_adh = self.params.demographics[agent_race][
                        agent_so
                    ].haart.adherence
                    if self.runRandom.random() < haart_adh:
                        adherence = 5
                    else:
                        adherence = self.runRandom.randint(1, 4)

                    # Add agent to HAART class set, update agent params
                    agent.haart = True
                    agent.intervention_ever = True
                    agent.haart_adherence = adherence
                    agent.haart_time = time
                    self.Trt_ART_agentSet.add_agent(agent)

            # Go off HAART
            elif (
                agent_haart
                and self.runRandom.random()
                < self.params.demographics[agent_race][agent_so].haart.discontinue
            ):
                agent.haart = False
                agent.haart_adherence = 0
                agent.haart_time = 0
                self.Trt_ART_agentSet.remove_agent(agent)

    def _discont_PrEP(self, agent: Agent, force: bool = False):
        # Agent must be on PrEP to discontinue PrEP
        assert agent.prep

        # If force flag set, auto kick off prep.
        if force:
            self.Trt_PrEP_agentSet.remove_agent(agent)
            self.prep_agents[agent.race][agent.so] -= 1
            agent.prep = False
            agent.prep_reason = []

        # else if agent is on PrEP, see if they should discontinue
        else:
            if (
                self.runRandom.random()
                < self.params.demographics[agent.race][agent.so].prep.discontinue
            ):
                self.Trt_PrEP_agentSet.remove_agent(agent)
                self.prep_agents[agent.race][agent.so] -= 1

                if "Oral" in self.params.prep.type:
                    agent.prep = False
                    agent.prep_reason = []
            else:  # if not discontinue, see if its time for a new shot. # REVIEWED what is this logic doing? This decrements, then update_prep_load increments - sarah to review with max
                if agent.prep_last_dose > 2:
                    agent.prep_last_dose = -1

        if self.params.prep.type == "Inj":
            agent.update_prep_load(self.params)

    def advance_vaccine(self, agent: Agent, time: int, vaxType: str):
        """
        :Purpose:
            Progress vaccine. Agents may receive injection or progress in time since injection.

        :Input:
            agent: Agent
            time: int

        :Output:
            none
        """
        if agent.vaccine:
            agent.vaccine_time += 1
            if (
                self.params.vaccine.booster
                and agent.vaccine_time
                == self.params.demographics[agent.race][
                    agent.so
                ].vaccine.booser.interval
                and self.runRandom.random()
                < self.params.demographics[agent.race][agent.so].vaccine.booster.prob
            ):
                agent.vaccinate(vaxType)

        elif time == self.params.vaccine.start:
            if (
                self.runRandom.random()
                < self.params.demographics[agent.race][agent.so].vaccine.prob
            ):
                agent.vaccinate(vaxType)

    def _initiate_PrEP(self, agent: Agent, time: int, force: bool = False):
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

        def _enrollPrEP(self, agent: Agent):
            agent.prep = True
            agent.intervention_ever = True
            self.Trt_PrEP_agentSet.add_agent(agent)
            self.newPrEPagents.add_agent(agent)

            self.prep_agents[agent.race][agent.so] += 1

            if (
                self.runRandom.random()
                < self.params.demographics[agent.race][agent.so].prep.adherence
            ):
                agent.prep_adherence = 1
            else:
                agent.prep_adherence = 0

            # set PrEP load and dosestep for PCK
            if "Inj" in self.params.prep.type:
                agent.prep_load = self.params.prep.peak_load
                agent.prep_last_dose = 0

        # agent must exist
        assert agent is not None

        # Prep only valid for agents not on prep and are HIV negative
        if agent.prep or agent.hiv:
            return

        # Determine probability of HIV treatment
        agent_race = agent.race

        if force:
            _enrollPrEP(self, agent)
        else:
            if self.params.prep.target_model == "Racial":
                numPrEP_agents = sum(self.prep_agents[agent_race].values())
            else:
                numPrEP_agents = self.Trt_PrEP_agentSet.num_members()

            if self.params.prep.target_model in ("Incar", "IncarHR"):
                if self.runRandom.random() < self.params.prep.target:
                    _enrollPrEP(self, agent)
                return None
            elif self.params.prep.target_model == "Racial":
                all_HIV_agents = set(self.All_agentSet.subset["HIV"].members)
                all_race = set(
                    self.All_agentSet.subset["Race"].subset[agent.race].members
                )
                HIV_agents = len(all_HIV_agents & all_race)
                # print("HIV agents", HIV_agents, "totHIV", len(all_HIV_agents))
                target_PrEP = (
                    int(
                        self.All_agentSet.subset["Race"]
                        .subset[agent.race]
                        .num_members()
                    )
                    - HIV_agents
                ) * self.params.demographics[agent.race][agent.so].prep.coverage

            else:
                target_PrEP = int(
                    (
                        self.All_agentSet.num_members()
                        - self.All_agentSet.subset["HIV"].num_members()
                    )
                    * self.params.prep.target
                )

            if self.params.prep.target_model in ("Incar", "IncarHR"):
                if self.runRandom.random() < self.params.prep.target:
                    _enrollPrEP(self, agent)
            elif (
                numPrEP_agents < target_PrEP
                and time >= self.params.prep.start
                and agent.PrEP_eligible(self.params.prep.target_model)
            ):
                _enrollPrEP(self, agent)

    def _progress_to_AIDS(self, agent: Agent):
        """
        :Purpose:
            Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        if not agent.hiv:
            raise ValueError("AIDS only valid for HIV agents!agent:%s" % str(agent.id))

        # REVIEWED Why do we check for not HAART, but then get HAART adherance? - Sarah to ask Max
        if not agent.haart:
            adherenceStat = agent.haart_adherence
            p = prob.adherence_prob(adherenceStat)

            if self.runRandom.random() < p * self.params.calibration.aids_progression:
                agent.aids = True
                self.HIV_AIDS_agentSet.add_agent(agent)

    def _die_and_replace(self):

        """
        :Purpose:
            Let agents die and replace the dead agent with a new agent randomly.
        """
        # die stage
        for agent in self.All_agentSet.members:

            # agent incarcerated, don't evaluate for death
            if agent.incar:
                continue

            # death rate per 1 person-month
            p = (
                prob.get_death_rate(
                    agent.hiv, agent.aids, agent.race, agent.haart_adherence
                )
                * self.params.calibration.mortality
            )

            if self.runRandom.random() < p:
                self.deathSet.append(agent)

                # End all existing relationships
                for rel in agent.relationships:
                    rel.progress(forceKill=True)

                    self.Relationships.remove(rel)

                # Remove agent node and edges from network graph
                self.get_Graph().remove_node(agent)

        # replace stage
        for agent in self.deathSet:
            # Remove agent from agent class and sub-sets
            self.All_agentSet.remove_agent(agent)

            new_agent = self.create_agent(agent.race, agent.so)
            self.add_agent_to_pop(new_agent)
            self.get_Graph().add_node(new_agent)
