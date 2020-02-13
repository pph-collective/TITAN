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


from .agent import Agent_set, Agent, Relationship
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
        self.NewInfections = Agent_set("NewInfections")
        self.NewDiagnosis = Agent_set("NewDiagnosis")
        self.NewIncarRelease = Agent_set("NewIncarRelease")
        self.NewHRrolls = Agent_set("NewHRrolls")
        self.newPrEPagents = Agent_set("NewPrEPagents")

        self.totalDiagnosis = 0
        self.treatmentEnrolled = False

        self.PrEPagents = {
            "BLACK": {"MSM": 0, "HF": 0, "HM": 0},
            "WHITE": {"MSM": 0, "HF": 0, "HM": 0},
        }
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

    def run(self):
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
                printer(run_id, t, self.runseed, self.popseed, self.netseed, stat)

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
            firstHIV = self.runRandom.choice(self.DU_IDU_agentSet._members)
            for i in range(numPartners):
                self.update_agent_partners(self.get_Graph(), firstHIV)
            self._become_HIV(firstHIV)

        run_id = uuid.uuid4()

        burnSimulation(params.model.burn_duration)

        print("\n === Begin Simulation Run ===")
        if params.outputs.draw_figures:
            nNodes = self.G.number_of_nodes()
            self.visualize_network(
                coloring=params.drawFigureColor,
                node_size=5000.0 / nNodes,
                curtime=0,
                iterations=10,
                label="Seed" + str(self.runseed),
            )

        if params.outputs.calc_component_stats:
            ao.print_components(
                run_id, 0, self.runseed, self.popseed, self.netseed, get_components()
            )

        print("\t===! Start Main Loop !===")

        # dictionary to hold results over time
        stats = {}

        # If we are using an agent zero method, create agent zero.
        if self.params.features.agent_zero:
            makeAgentZero(4)

        if self.params.outputs.edge_list:
            path = "results/network/Edgelist_t0.txt"
            self.write_G_edgelist(path)

        for t in range(1, self.tmax + 1):
            print(f"\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t.: TIME {t}")
            if (
                params.outputs.draw_figures
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
                self.params
            )
            print_stats(stats[t], run_id)

            print(("Number of relationships: %d" % len(self.Relationships)))
            self.All_agentSet.print_subsets()

            self.totalDiagnosis += len(self.NewDiagnosis._members)
            if (
                self.totalDiagnosis > self.params.needle_exchange.init_at_pop
                and not self.treatmentEnrolled
            ):
                self._enroll_treatment()

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
                    )
                if self.params.outputs.edge_list:
                    path = f"results/network/Edgelist_t{t}.txt"
                    self.write_G_edgelist(path)

        return stats

    def _update_AllAgents(self, time: int, burn: bool = False):
        """
        :Purpose:
            Update IDU agents:
            For each agent:
                1 - determine agent type
                2 - get partners
                3 - agent interacts with partners
                5 - VCT (Voluntsry Counseling and Testing)
                6 - if IDU: SEP, treatment
                7 - if HIV: HAART, AIDS

        :Input:
            agent, time

        :Output:
            none
        """
        if time > 0 and self.params.features.static_n is False:
            self.update_partner_assignments(self.get_Graph())

        for rel in self.Relationships:
            # If in burn, ignore interactions
            if not burn:
                self.pop_randomract(rel)

            # If static network, ignore relationship progression
            if not self.params.features.static_n:
                if rel.progress():
                    g = self.get_Graph()
                    if g.has_edge(rel._ID1, rel._ID2):
                        g.remove_edge(rel._ID1, rel._ID2)

                    self.Relationships.remove(rel)
                    del rel

        if self.params.features.high_risk:
            for agent in self.highrisk_agentsSet.iter_agents():
                if agent._highrisk_time > 0:
                    agent._highrisk_time -= 1
                    if (
                        agent._SO == "HM"
                        and self.params.features.prep
                        and (self.params.prep.target_model in ("HR", "IncarHR"))
                    ):
                        for part in agent._partners:
                            if not (part._HIV_bool or part.vaccine_bool):
                                self._initiate_PrEP(part, time)
                else:
                    self.highrisk_agentsSet.remove_agent(agent)
                    agent._highrisk_bool = False

                    if (
                        self.params.features.incar
                    ):  # TO_REVIEW why does this check hm/hf and then subtract things - could this be more generic?
                        if agent._SO == "HM":
                            agent._mean_num_partners -= (
                                self.params.high_risk.partner_scale
                            )
                        elif agent._SO == "HF":
                            agent._mean_num_partners -= (
                                self.params.high_risk.partner_scale
                            )

        for agent in self.All_agentSet.iter_agents():
            agent._timeAlive += 1

            if self.params.features.incar:
                self._incarcerate(agent, time)

            if agent._MSMW and self.runRandom.random() < params.HIV_MSMW:
                self._become_HIV(agent)

            if agent._HIV_bool:
                # If in burnin, ignore HIV
                if burn:
                    if (
                        agent._incar_treatment_time >= 1
                    ):  # TO_REVIEW why does this decrement only during burn?
                        agent._incar_treatment_time -= 1

                else:
                    self._HIVtest(agent, time)
                    self._progress_to_AIDS(agent)

                    if self.params.features.haart:
                        self._update_HAART(agent, time)
                        agent._HIV_time += 1
            else:
                if self.params.features.prep:
                    if time >= self.params.prep.start:
                        if agent._PrEP_bool:
                            self._discont_PrEP(agent)
                        elif self.params.prep.target_model == "RandomTrial":
                            pass
                        elif agent.PrEP_eligible():
                            self._initiate_PrEP(agent, time)
                    if (
                        self.params.features.vaccine
                        and not agent._PrEP_bool
                        and not burn
                    ):
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
                            if (not ag._HIV_bool) and (not ag._PrEP_bool):
                                ag._treatment_bool = True
                                if (
                                    self.runRandom.random() < self.params.prep.target
                                    and not agent.vaccine_bool
                                ):
                                    self._initiate_PrEP(ag, time, force=True)
                print(("Total agents in trial: ", totNods))

    def _agents_interact(self, rel: Relationship) -> bool:
        """
        :Purpose:
            Let IDU agent interact with a partner.
            Update IDU agents:
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
        if rel._ID1._incar_bool or rel._ID2._incar_bool:
            return False

        if (
            rel._ID1._HIV_bool and not rel._ID2._HIV_bool
        ):  # Agent 1 is HIV, partner is succept
            agent = rel._ID1
            partner = rel._ID2
        elif (
            not rel._ID1._HIV_bool and rel._ID2._HIV_bool
        ):  # If agent_2 is HIV agen1 is not, agent_2 is HIV, agent_1 is succept
            agent = rel._ID2
            partner = rel._ID1
        else:  # neither agent is HIV or both are
            return False

        rel_sex_possible = sex_possible(agent._SO, partner._SO, self.params)
        partner_drug_type = partner._DU
        agent_drug_type = agent._DU

        if partner_drug_type == "IDU" and agent_drug_type == "IDU":
            # Injection is possible
            # If agent is on post incar HR treatment to prevent IDU behavior, pass IUD infections

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

        elif partner_drug_type in ["NIDU", "NDU"] or agent_drug_type in ["NIDU", "NDU"]:
            if rel_sex_possible:
                self._sex_transmission(rel)
            else:
                return False
        else:  # REVIEWED - sanity test, with params re-write this logic/check can move there
            raise ValueError("Agents must be either IDU, NIDU, or ND")

        return True  # if didn't short circuit, agents interacted

    def _needle_transmission(self, agent: Agent, partner: Agent):
        """
        :Purpose:
            Simulate random transmission of HIV between two IDU agents
            through needle.\n
            Agent must by HIV+ and partner not.

        :Input:
            agents : int
            partner : int
        :Output: -
        """

        assert agent._HIV_bool
        assert not partner._HIV_bool
        assert agent._DU == "IDU"
        assert partner._DU == "IDU"

        agent_race = agent._race
        agent_sex_type = agent._SO

        # REVIEWED why is the mean number of sex acts for a class multiplied by needle calibration? - change to num_needle_acts
        MEAN_N_ACTS = (
            self.params.demographics[agent_race][agent_sex_type].num_needle_acts
            * self.params.calibration.needle.act
        )
        share_acts = round(poisson.rvs(MEAN_N_ACTS, size=1)[0])

        if agent._SNE_bool:  # safe needle exchange - minimal sharing
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

        if rel._ID1._HIV_bool:
            agent = rel._ID1
            partner = rel._ID2
        elif rel._ID2._HIV_bool:
            agent = rel._ID2
            partner = rel._ID1
        else:
            raise ValueError("rel must have an agent with HIV")

        # HIV status of agent and partner
        # Everything from here is only run if one of them is HIV+
        if partner._HIV_bool:
            return

        # unprotected sex probabilities for primary partnerships
        MSexActs = (
            agent.get_number_of_sexActs(self.runRandom, self.params)
            * self.params.calibration.sex.act
        )
        T_sex_acts = round(poisson.rvs(MSexActs, size=1)[0])

        # Get condom usage
        if self.params.high_risk.condom_use_type == "Race":
            p_SafeSex = self.params.demographics[agent._race][agent._SO].safe_sex
        else:
            p_SafeSex = prob.safe_sex(rel._total_sex_acts)

        # Reduction of risk acts between partners for condom usage
        U_sex_acts = T_sex_acts
        for n in range(U_sex_acts):
            if self.runRandom.random() < p_SafeSex:
                U_sex_acts -= 1

        if U_sex_acts >= 1:
            # agent is HIV+
            rel._total_sex_acts += U_sex_acts
            ppAct = agent.get_transmission_probability("SEX", self.params)

            # Reduction of transmissibility for acts between partners for PrEP adherence
            if agent._PrEP_bool or partner._PrEP_bool:
                if "Oral" in self.params.prep.type:  # params.prep.type == "Oral":
                    if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                        ppAct = ppAct * (
                            1.0 - self.params.prep.efficacy.adherent
                        )  # 0.04
                    else:
                        ppAct = ppAct * (
                            1.0 - self.params.prep.efficacy.non_adherant
                        )  # 0.24

                elif "Inj" in self.params.prep.type:
                    ppActReduction = (
                        -1.0 * np.exp(-5.528636721 * partner._PrEP_load) + 1
                    )
                    if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                        ppAct = ppAct * (1.0 - ppActReduction)  # 0.04

            if partner.vaccine_bool:
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
        if not agent._HIV_bool:
            agent._HIV_bool = True
            agent._HIV_time = 1
            agent.vaccine_bool = False
            self.NewInfections.add_agent(agent)
            self.HIV_agentSet.add_agent(agent)

        if agent._PrEP_bool:
            self._discont_PrEP(agent, force=True)

    def _enroll_treatment(
        self
    ):  # TO_REVIEW can this be changed to enroll_needle_exchange?
        """
        :Purpose:
            Enroll IDU agents in needle exchange
        """
        print(("\n\n!!!!Engaginge treatment process"))
        self.treatmentEnrolled = True
        for agent in self.All_agentSet.iter_agents():
            if (
                self.runRandom.random() < self.params.needle_exchange.coverage
                and agent._DU == "IDU"
            ):
                agent._SNE_bool = True

    # REVIEWED this isn't used anywhere, but should be! _incarcerate makes things high risk and should reference this
    def _becomeHighRisk(self, agent: Agent, duration: int = None):

        if agent not in self.highrisk_agentsSet._members:
            self.highrisk_agentsSet.add_agent(agent)

        if not agent._everhighrisk_bool:
            self.NewHRrolls.add_agent(agent)

        agent._highrisk_bool = True
        agent._everhighrisk_bool = True

        if duration is not None:
            agent._highrisk_time = duration
        else:
            agent._highrisk_time = params.HR_M_dur

    def _incarcerate(self, agent: Agent, time: int):
        """
        :Purpose:
            To incarcerate an agent or update their incarceration variables

        :Input:
            agent : int
            time : int

        """
        hiv_bool = agent._HIV_bool

        if agent._incar_bool:
            agent._incar_time -= 1

            if agent._incar_time == 0:  # FREE AGENT
                self.incarcerated_agentSet.remove_agent(agent)
                self.NewIncarRelease.add_agent(agent)
                agent._incar_bool = False
                agent._ever_incar_bool = True
                if (
                    self.params.features.incar
                ):  # TO_REVIEW this was looking for model = incar which is broader - does this substitution make sense?
                    if (
                        not agent._highrisk_bool and self.params.features.high_risk
                    ):  # If behavioral treatment on and agent HIV, ignore HR period.
                        if (
                            self.params.incar.treatment.high_risk
                            and hiv_bool
                            and (time >= self.params.incar.treatment.start)
                        ):
                            pass
                        else:  # Else, become high risk
                            self.highrisk_agentsSet.add_agent(agent)
                            if not agent._everhighrisk_bool:
                                self.NewHRrolls.add_agent(agent)

                            agent._mean_num_partners = (
                                agent._mean_num_partners
                                + self.params.high_risk.partner_scale
                            )
                            agent._highrisk_bool = True
                            agent._everhighrisk_bool = True
                            agent._highrisk_time = self.params.HR_M_dur

                    if (
                        self.params.incar.treatment.ric
                        or self.params.incar.treatment.high_risk
                        or self.params.incar.treatment.idu
                    ) and (time >= self.params.incar.treatment.start):
                        agent._incar_treatment_time = (
                            self.params.incar.treatment.duration
                        )

                    if hiv_bool:
                        if agent._HAART_bool:
                            if (
                                self.runRandom.random()
                                > self.params.incar.haart.discontinue
                            ):  # 12% remain surpressed
                                pass

                            else:
                                agent._HAART_bool = False
                                agent._HAART_adh = 0
                                self.Trt_ART_agentSet.remove_agent(agent)

                            # END FORCE

        elif self.runRandom.random() < (
            self.params.demographics[agent._race][agent._SO].incar.prob
            * (1 + (hiv_bool * 4))
            * self.params.calibration.incarceration
        ):
            # REVIEWED what about other sex types? -needs to be generalized - Sarah meeting with someone
            if agent._SO == "HF":
                jailDuration = prob.HF_jail_duration
            elif agent._SO == "HM":
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
                if not agent._tested:
                    if self.runRandom.random() < self.params.incar.hiv.dx:
                        agent._tested = True
                else:  # Then tested and HIV, check to enroll in ART
                    if self.runRandom.random() < self.params.incar.haart.prob:
                        tmp_rnd = self.runRandom.random()
                        haart_adh = self.params.incar.haart.adherence
                        if tmp_rnd < haart_adh:
                            adherence = 5
                        else:
                            adherence = self.runRandom.randint(1, 4)

                        # Add agent to HAART class set, update agent params
                        agent._HAART_bool = True
                        agent._HAART_adh = adherence
                        agent._HAART_time = time
                        self.Trt_ART_agentSet.add_agent(agent)

            agent._incar_bool = True
            agent._incar_time = timestay
            self.incarcerated_agentSet.add_agent(agent)

            # PUT PARTNERS IN HIGH RISK
            for partner in agent._partners:
                if not partner._highrisk_bool:
                    if self.runRandom.random() < self.params.high_risk.proportion:
                        if not partner._highrisk_bool:
                            self.highrisk_agentsSet.add_agent(partner)
                            if not partner._everhighrisk_bool:
                                self.NewHRrolls.add_agent(partner)
                            partner._mean_num_partners += (
                                self.params.high_risk.partner_scale
                            )  # 32.5 #2 + 3.25 from incar HR
                            partner._highrisk_bool = True
                            partner._everhighrisk_bool = True
                            partner._highrisk_time = self.params.high_risk.sex_based[partner._SO].duration
                if self.params.features.prep and (
                    self.params.prep.target_model in ("Incar", "IncarHR")
                ):
                    # Atempt to put partner on prep if less than probability
                    if not partner._HIV_bool and not agent.vaccine_bool:
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
        sex_type = agent._SO
        race_type = agent._race
        tested = agent._tested

        def diagnose(agent):
            agent._tested = True
            self.NewDiagnosis.add_agent(agent)
            self.Trt_Tstd_agentSet.add_agent(agent)
            if (
                self.params.features.partner_tracing
            ):  # TODO fix this logic; should get partnerTraced and then lose it after
                # For each partner, determine if found by partner testing # TO_REVIEW how to deal with setting based logic - moved the intended functionality to a flagged feature
                for ptnr in agent._partners:
                    if ptnr._HIV_bool and not ptnr._tested:
                        ptnr.partnerTraced = True
                        ptnr.traceTime = time + 1

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
                agent.partnerTraced
                and self.runRandom.random() < 0.87
                and agent.traceTime == time
            ):
                diagnose(agent)

        agent.partnerTraced = False

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
        assert agent._HIV_bool

        agent_haart = agent._HAART_bool
        agent_race = agent._race
        agent_so = agent._SO

        # Determine probability of HIV treatment
        if time >= 0 and agent._tested:
            # Go on HAART
            if not agent_haart and agent._HAART_time == 0:
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
                    agent._HAART_bool = True
                    agent._HAART_adh = adherence
                    agent._HAART_time = time
                    self.Trt_ART_agentSet.add_agent(agent)

            # Go off HAART
            elif (
                agent_haart
                and self.runRandom.random()
                < self.params.demographics[agent_race][agent_so].haart.discontinue
            ):
                if not (
                    agent._incar_treatment_time > 0 and self.params.incar.treatment.ric
                ):
                    agent._HAART_bool = False
                    agent._HAART_adh = 0
                    agent._HAART_time = 0
                    self.Trt_ART_agentSet.remove_agent(agent)

    def _discont_PrEP(self, agent: Agent, force: bool = False):
        # Agent must be on PrEP to discontinue PrEP
        assert agent._PrEP_bool

        # If force flag set, auto kick off prep.
        if force:
            self.Trt_PrEP_agentSet.remove_agent(agent)
            self.PrEPagents[agent._race][agent._SO] -= 1
            agent._PrEP_bool = False
            agent._PrEP_reason = []

        # else if agent is on PrEP, see if they should discontinue
        else:
            if (
                self.runRandom.random()
                < self.params.demographics[agent._race][agent._SO].prep.discontinue
            ):
                self.Trt_PrEP_agentSet.remove_agent(agent)
                self.PrEPagents[agent._race][agent._SO] -= 1

                if "Oral" in self.params.prep.type:
                    agent._PrEP_bool = False
                    agent._PrEP_reason = []
            else:  # if not discontinue, see if its time for a new shot. # REVIEWED what is this logic doing? This decrements, then update_PrEP_load increments - sarah to review with max
                if agent._PrEP_lastDose > 2:
                    agent._PrEP_lastDose = -1

        if self.params.prep.type == "Inj":
            agent.update_PrEP_load()

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
        if agent.vaccine_bool:
            agent.vaccine_time += 1
            if (
                self.params.vaccine.booster
                and agent.vaccine_time
                == self.params.demographics[agent._race][
                    agent._SO
                ].vaccine.booser.interval
                and self.runRandom.random()
                < self.params.demographics[agent._race][agent._SO].vaccine.booster.prob
            ):
                agent.vaccinate(vaxType)

        elif time == self.params.vaccine.start:
            if (
                self.runRandom.random()
                < self.params.demographics[agent._race][agent._SO].vaccine.prob
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
            agent._PrEP_bool = True
            self.Trt_PrEP_agentSet.add_agent(agent)
            self.newPrEPagents.add_agent(agent)

            self.PrEPagents[agent._race][agent._SO] += 1

            tmp_rnd = self.runRandom.random()
            # if params.setting == "AtlantaMSM":  # TO_REVIEW setting based logic
            #     if (
            #         tmp_rnd
            #         < self.params.demographics[agent._race][agent._SO].prep.adherence
            #     ):
            #         agent._PrEP_adh = 1
            #     else:
            #         agent._PrEP_adh = 0
            # else:
            if tmp_rnd < self.params.prep.adherence:
                agent._PrEP_adh = 1
            else:
                agent._PrEP_adh = 0

            # set PrEP load and dosestep for PCK
            if "Inj" in self.params.prep.type:
                agent._PrEP_load = self.params.prep.peak_load
                agent._PrEP_lastDose = 0

        # agent must exist
        assert agent is not None

        # Prep only valid for agents not on prep and are HIV negative
        if agent._PrEP_bool or agent._HIV_bool:
            return

        # Determine probability of HIV treatment
        agent_race = agent._race

        if force:
            _enrollPrEP(self, agent)
        else:
            if self.params.prep.target_model == "Racial":
                numPrEP_agents = sum(self.PrEPagents[agent_race].values())
            else:
                numPrEP_agents = self.Trt_PrEP_agentSet.num_members()

            if self.params.prep.target_model in ("Incar", "IncarHR"):
                if self.runRandom.random() < self.params.prep.target:
                    _enrollPrEP(self, agent)
                return None
            elif self.params.prep.target_model == "Racial":
                all_HIV_agents = set(self.All_agentSet._subset["HIV"]._members)
                all_race = set(
                    self.All_agentSet._subset["Race"]._subset[agent._race]._members
                )
                HIV_agents = len(all_HIV_agents & all_race)
                # print("HIV agents", HIV_agents, "totHIV", len(all_HIV_agents))
                target_PrEP = (
                    int(
                        self.All_agentSet._subset["Race"]
                        ._subset[agent._race]
                        .num_members()
                    )
                    - HIV_agents
                ) * self.params.demographics[agent._race][agent._SO].prep.coverage

            else:
                target_PrEP = int(
                    (
                        self.All_agentSet.num_members()
                        - self.All_agentSet._subset["HIV"].num_members()
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
        if not agent._HIV_bool:
            raise ValueError("AIDS only valid for HIV agents!agent:%s" % str(agent._ID))

        # REVIEWED Why do we check for not HAART, but then get HAART adherance? - Sarah to ask Max
        if not agent._HAART_bool:
            adherenceStat = agent._HAART_adh
            p = prob.adherence_prob(adherenceStat)

            if self.runRandom.random() < p * self.params.calibration.aids_progression:
                agent._AIDS_bool = True
                self.HIV_AIDS_agentSet.add_agent(agent)

    def _die_and_replace(self):

        """
        :Purpose:
            Let agents die and replace the dead agent with a new agent randomly.
        """
        # die stage
        for agent in self.All_agentSet._members:

            # agent incarcerated, don't evaluate for death
            if agent._incar_bool:
                continue

            # death rate per 1 person-month
            p = (
                prob.get_death_rate(
                    agent._HIV_bool, agent._AIDS_bool, agent._race, agent._HAART_adh
                )
                * self.params.calibration.mortality
            )

            if self.runRandom.random() < p:
                self.deathSet.append(agent)

                # End all existing relationships
                for rel in agent._relationships:
                    rel.progress(forceKill=True)

                    self.Relationships.remove(rel)

                # Remove agent node and edges from network graph
                self.get_Graph().remove_node(agent)

        # replace stage
        for agent in self.deathSet:
            # Remove agent from agent class and sub-sets
            self.All_agentSet.remove_agent(agent)

            new_agent = self.create_agent(agent._race, agent._SO)
            self.add_agent_to_pop(new_agent)
            self.get_Graph().add_node(new_agent)
