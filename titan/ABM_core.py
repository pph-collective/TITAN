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


from .agent import Agent_set, Agent, Relationship
from .network_graph_tools import NetworkClass
from . import analysis_output as ao
from . import probabilities as prob
from . import params  # type: ignore
from .ABM_partnering import sex_possible


class HIVModel(NetworkClass):
    """
    :Purpose:
        This is the core class used to simulate
        the spread of HIV and drug use in one MSA
        (Metropolitan Statistical Area).

    :Input:
        N : int
            Number of agents. Default: 1000
        tmax: int
            Number of simulation steps (years).
        runseed: int
            random seed for running the model
        popseed: int
            random seed for initalizing population
        netseed: int
            random seed for initializing network
        network_type: str
            type of network (e.g. "scale_free")
    """

    def __repr__(self):
        returnStr = "\n"
        returnStr += "Seed: %d\n" % (self.runseed)
        returnStr += "Npop: %d\n" % (params.N_POP)
        returnStr += "Time: %d\n" % (params.TIME_RANGE)
        returnStr += "Mode: %s\n" % (params.model)

        return returnStr

    def __init__(
        self,
        N: int,
        tmax: int,
        runseed: int,
        popseed: int,
        netseed: int,
        network_type: str = "",
    ):
        # Ensure param variables are defined. For backwards compatibility with params.py files
        bc_attrs = [
            "drawEdgeList",
            "inc_treat_HRsex_beh",
            "inc_treat_IDU_beh",
            "calcNetworkStats",
        ]
        for attr in bc_attrs:
            if not hasattr(params, attr):
                setattr(params, attr, False)

        if type(tmax) is not int:
            raise ValueError("Number of time steps must be integer")
        else:
            self.tmax = tmax

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

        self.runseed = get_check_rand_int(runseed)
        self.popseed = get_check_rand_int(popseed)
        self.netseed = get_check_rand_int(netseed)

        print("=== Begin Initialization Protocol ===\n")

        NetworkClass.__init__(
            self,
            N=N,
            network_type=network_type,
            popSeed=self.popseed,
            netSeed=self.netseed,
        )

        print("\n\tCreating lists")
        # Other lists / dictionaries

        self.NewInfections = Agent_set("NewInfections")
        self.NewDiagnosis = Agent_set("NewDiagnosis")
        self.NewIncarRelease = Agent_set(
            "NewIncarRelease"
        )  # REVIEW: ask wgoedel if this is needed (doesn't print)
        self.NewHRrolls = Agent_set("NewHRrolls")  # REVIEW: ditto above

        self.totalDiagnosis = 0
        self.treatmentEnrolled = False

        self.newPrEPagents = Agent_set("NewPrEPagents")  # REVIEW is this ever used
        self.PrEPagents = {
            "BLACK": {"MSM": 0, "HF": 0, "HM": 0},
            "WHITE": {"MSM": 0, "HF": 0, "HM": 0},
        }
        # Set seed format. 0: pure random, -1: Stepwise from 1 to nRuns, else: fixed value
        print(("\tRun seed was set to:", runseed))
        self.runRandom = Random(runseed)
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
            for report in params.reports:
                printer = getattr(ao, report)
                printer(run_id, t, self.runseed, self.popseed, self.netseed, stat)

        def get_components():
            return list(
                self.G.subgraph(c).copy()
                for c in nx.connected_components(self.G)
                if len(c) >= params.minComponentSize
            )

        def burnSimulation(burnDuration: int):
            print(
                ("\n === Burn Initiated for {} timesteps ===".format(burnDuration + 1))
            )
            for t in range(0, burnDuration + 1):
                self._update_AllAgents(t, burn=True)
                if params.flag_DandR:
                    self._die_and_replace(t)
            self.All_agentSet.print_subsets()
            print(("\tBurn Cuml Inc:\t{}".format(self.NewInfections.num_members())))
            self.NewInfections.clear_set()
            self.NewDiagnosis.clear_set()
            self.NewHRrolls.clear_set()
            self.NewIncarRelease.clear_set()
            self.newPrEPagents.clear_set()

            self.deathSet = []
            print(" === Simulation Burn Complete ===")

        run_id = uuid.uuid4()

        burnSimulation(params.burnDuration)

        print("\n === Begin Simulation Run ===")
        if params.drawFigures:  # REVIEW: is this ever used? seems deprecated
            nNodes = self.G.number_of_nodes()
            self.visualize_network(
                coloring=params.drawFigureColor,
                node_size=5000.0 / nNodes,
                curtime=0,
                iterations=10,
                label="Seed" + str(self.runseed),
            )

        if params.calcComponentStats:
            ao.print_components(
                run_id, 0, self.runseed, self.popseed, self.netseed, get_components()
            )

        def makeAgentZero(numPartners: int):
            firstHIV = self.runRandom.choice(self.DU_IDU_agentSet._members)
            for i in range(numPartners):
                self.update_agent_partners(self.get_Graph(), firstHIV)
            self._become_HIV(firstHIV, 0)

        print("\t===! Start Main Loop !===")

        # dictionary to hold results over time
        stats = {}

        # If we are using an agent zero method, create agent zero.
        if params.flag_agentZero:
            makeAgentZero(4)

        if params.drawEdgeList:
            print("Drawing network edge list to file")
            fh = open("results/network/Edgelist_t{}.txt".format(0), "wb")
            self.write_G_edgelist(fh)
            fh.close()

        for t in range(1, self.tmax + 1):
            print(f"\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t.: TIME {t}")
            if params.drawFigures and t % params.intermPrintFreq == 0:
                self.visualize_network(
                    coloring=params.drawFigureColor,
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

            if params.flag_DandR:
                self._die_and_replace(t)

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
                self.LAI_agentSet,
                self.oralPrEP_agentSet,
            )
            print_stats(stats[t], run_id)

            print(("Number of relationships: %d" % len(self.Relationships)))
            self.All_agentSet.print_subsets()

            self.totalDiagnosis += len(self.NewDiagnosis._members)
            if (
                self.totalDiagnosis > params.initTreatment
                and not self.treatmentEnrolled
            ):
                self._enroll_treatment(t)

            # RESET counters for the next time step
            self.deathSet = []
            self.NewInfections.clear_set()
            self.NewDiagnosis.clear_set()
            self.NewHRrolls.clear_set()
            self.NewIncarRelease.clear_set()
            self.newPrEPagents.clear_set()

            print((t % params.intermPrintFreq))  # REVIEW: is this broken?
            if t % params.intermPrintFreq == 0:
                if params.calcNetworkStats:
                    self.write_network_stats(t=t)
                if params.calcComponentStats:
                    ao.print_components(
                        run_id,
                        t,
                        self.runseed,
                        self.popseed,
                        self.netseed,
                        get_components(),
                    )
                if params.drawEdgeList:  # REVIEW: is this broken?
                    print("Drawing network edge list to file")
                    fh = open("results/network/Edgelist_t{}.txt".format(t), "wb")
                    self.write_G_edgelist(fh)
                    fh.close()

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
        if time > 0 and params.flag_staticN is False:
            self.update_partner_assignments(params.PARTNERTURNOVER, self.get_Graph())

        for rel in self.Relationships:
            # If in burn, ignore interactions
            if burn:
                pass
            else:
                self._agents_interact(time, rel)

            # If static network, ignore relationship progression
            if params.flag_staticN:
                pass
            else:
                if rel.progress():
                    g = self.get_Graph()
                    if g.has_edge(rel._ID1, rel._ID2):
                        g.remove_edge(rel._ID1, rel._ID2)

                    self.Relationships.remove(rel)
                    del rel

        if params.flag_high_risk:  # TODO: abstract this
            for tmpA in self.highrisk_agentsSet.iter_agents():
                if tmpA._highrisk_time > 0:
                    tmpA._highrisk_time -= 1
                    if (
                        tmpA._SO == "HM"
                        and params.flag_PrEP
                        and (
                            params.PrEP_target_model == "high_risk"
                            or params.PrEP_target_model == "incarcerated_high_risk"
                        )
                    ):
                        for part in tmpA._partners:
                            if not (part._HIV_bool or part.vaccine_bool):
                                self._initiate_PrEP(part, time)
                else:
                    self.highrisk_agentsSet.remove_agent(tmpA)
                    tmpA._highrisk_bool = False

                    if params.model == "Incar":  # TODO abstract this
                        if tmpA._SO == "HM":
                            tmpA._mean_num_partners -= params.HR_partnerScale
                        elif tmpA._SO == "HF":
                            tmpA._mean_num_partners -= params.HR_partnerScale

        for agent in self.All_agentSet.iter_agents():
            agent._timeAlive += 1
            if (
                params.flag_PCA
                and self.runRandom.random() < params.awarenessProb
                and not burn
            ):
                agent.awareness = True
                self.aware_agentSet.add_agent(agent)
                if self.runRandom.random() < params.PCA_PrEP:
                    self._initiate_PrEP(agent, time, force=True)
            if params.flag_incar:  # and not burn:
                self._incarcerate(agent, time)
            if agent._MSMW and self.runRandom.random() < params.HIV_MSMW:
                self._become_HIV(agent, 0)
            if agent._HIV_bool:
                # If in burnin, ignore HIV
                if burn:
                    if agent._incar_treatment_time >= 1:
                        agent._incar_treatment_time -= 1

                else:
                    self._HIVtest(agent, time)
                    self._progress_to_AIDS(agent)

                    if params.flag_ART:
                        self._update_HAART(agent, time)
                        agent._HIV_time += 1
            else:
                if params.flag_PrEP:
                    if time >= params.PrEP_startT:
                        if agent._PrEP_bool:
                            self._discont_PrEP(agent, time)
                        elif params.PrEP_target_model == "Clinical":
                            pass
                        elif params.PrEP_target_model == "RandomTrial":
                            pass
                        elif agent.PrEP_eligible(time) and not agent._PrEP_bool:
                            self._initiate_PrEP(agent, time)
                    if "Vaccine" in params.PrEP_type and not agent._PrEP_bool:
                        self.advance_vaccine(
                            agent, time, vaxType=params.vaccine_type, burn=burn
                        )
        if params.flag_PrEP and time >= params.PrEP_startT:
            if params.PrEP_target_model == "Clinical":  # REVIEW is this used?
                if time > params.PrEP_startT:
                    numPrEP_agents = self.Trt_PrEP_agentSet.num_members()
                    target_PrEP = int(
                        (
                            self.All_agentSet.num_members()
                            - self.All_agentSet._subset["HIV"].num_members()
                        )
                        * params.PrEP_Target
                    )
                    eligiblePool = [
                        ag
                        for ag in self.All_agentSet._subset["SO"]
                        ._subset["MSM"]
                        ._members
                        if (ag._PrEP_bool is False and ag._HIV_bool is False)
                    ]

                    while numPrEP_agents < target_PrEP:
                        numPrEP_agents = self.Trt_PrEP_agentSet.num_members()
                        target_PrEP = int(
                            (
                                self.All_agentSet.num_members()
                                - self.All_agentSet._subset["HIV"].num_members()
                            )
                            * params.PrEP_Target
                        )
                        selectedAgent = self._get_clinic_agent(
                            params.PrEP_clinic_cat, eligiblePool
                        )
                        if selectedAgent is not None:
                            eligiblePool.remove(
                                selectedAgent
                            )  # shouldn't be selected again
                            self._initiate_PrEP(selectedAgent, time)
            elif (
                "RandomTrial" in params.PrEP_target_model and time == params.PrEP_startT
            ):
                print("Starting random trial")
                components = list(
                    self.G.subgraph(c).copy()
                    for c in nx.connected_components(self.G)
                    if len(c) >= params.minComponentSize
                )
                totNods = 0
                print(
                    "Number of components",
                    len(
                        [
                            1
                            for comp in components
                            if comp.number_of_nodes() >= params.minComponentSize
                        ]
                    ),
                )
                for comp in components:
                    totNods += comp.number_of_nodes()

                    if self.runRandom.random() < 0.5:
                        # Component selected as treatment pod!
                        if not params.flag_PCA:
                            for ag in comp.nodes():
                                if (
                                    (ag._HIV_bool is False)
                                    and (ag._PrEP_bool is False)
                                    and not params.flag_PCA
                                ):
                                    ag._treatment_bool = True
                                    if (
                                        self.runRandom.random() < params.PrEP_Target
                                        and not ag.vaccine_bool
                                    ):
                                        self._initiate_PrEP(ag, time, force=True)
                        if params.pcaChoice == "eigenvector":
                            centrality = nx.eigenvector_centrality_numpy(self.G)
                            assert len(centrality) >= 1, "Empty centrality"
                            orderedCentrality = sorted(centrality, key=centrality.get)
                            intervention_agent = False
                            for ag in orderedCentrality:
                                if not ag._HIV_bool:
                                    self.PCA_agentSet.add_agent(ag)
                                    ag.awareness = True
                                    ag._PCA = 1
                                    self.aware_agentSet.add_agent(ag)
                                    intervention_agent = True
                                    break
                                else:
                                    pass
                            if not intervention_agent:
                                ag = orderedCentrality[0]
                                ag._PCA = -1

                        elif params.pcaChoice == "bridge":
                            all_bridges = list(
                                nx.bridges(comp)
                            )  # get a list of bridges
                            comp_agents = [
                                ag for ag, j in all_bridges if not ag._HIV_bool
                            ]  # all agents in bridge (first half)
                            if comp_agents:
                                chosen_agent = random.choice(
                                    comp_agents
                                )  # select change agent
                                self.aware_agentSet.add_agent(
                                    chosen_agent
                                )  # add to aware agents
                                chosen_agent.awareness = True  # make aware
                                self.PCA_agentSet.add_agent(
                                    chosen_agent
                                )  # add to PCA agents
                                chosen_agent._PCA = 1
                            elif all_bridges:
                                chosen_bridge = random.choice(
                                    list(all_bridges)
                                )  # not true change agent, just mark component
                                chosen_agent = random.choice(list(chosen_bridge))
                                chosen_agent._PCA = -1  # make change agent
                            else:
                                chosen_agent = random.choice(list(comp.nodes))
                                chosen_agent._PCA = -1

                        elif params.pcaChoice == "random":
                            chosen_agent = None
                            for ag in comp:
                                if not ag._HIV_bool:
                                    chosen_agent = ag
                                    chosen_agent._PCA = 1
                                    self.aware_agentSet.add_agent(
                                        chosen_agent
                                    )  # add to aware agents
                                    chosen_agent.awareness = True  # make aware
                                    self.PCA_agentSet.add_agent(
                                        chosen_agent
                                    )  # add to PCA agents
                                    break
                                else:
                                    pass
                            if not chosen_agent:
                                chosen_agent = random.choice(list(comp.nodes))
                                chosen_agent._PCA = -1

                print(("Total agents in trial: ", totNods))
                print("Number of Change Agents:", self.PCA_agentSet.num_members())

    def _agents_interact(self, time: int, rel: Relationship):
        """f
            Let IDU agent interact with a partner.
            Update IDU agents:
                1 - determine transition type
                2 - Injection rules
                3 - Sex rules
                4 - HIV transmission

        :Input:

            time : int

            rel : Relationship

            rand_gen : random number generator

        Output:
            none

        """
        # If either agent is incarcerated, skip their interaction
        if rel._ID1._incar_bool or rel._ID2._incar_bool:
            return

        if params.flag_PCA:
            if rel._rel_type == "sexOnly" and rel._duration != rel._total_duration:
                pass
            else:
                self._pca_interaction(rel, time)

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
            return

        rel_sex_possible = sex_possible(agent._SO, partner._SO)
        partner_drug_type = partner._DU
        agent_drug_type = agent._DU

        if partner_drug_type == "IDU" and agent_drug_type == "IDU":
            # Injection is possible
            # If agent is on post incar HR treatment to prevent IDU behavior, pass IUD infections
            if agent._incar_treatment_time > 0 and params.inc_treat_IDU_beh:
                return

            elif rel_sex_possible:
                # Sex is possible
                rv = self.runRandom.random()  # REVIEW after bond types established
                if rv < 0.25:  # Needle only (60%)
                    self._needle_transmission(agent, partner, time)
                else:  # Both sex and needle (20%)
                    self._needle_transmission(agent, partner, time)
                    self._sex_transmission(time, rel)
            else:
                # Sex not possible, needle only
                self._needle_transmission(agent, partner, time)

        elif partner_drug_type in ["NIDU", "NDU"] or agent_drug_type in ["NIDU", "NDU"]:
            if rel_sex_possible:
                self._sex_transmission(time, rel)
            else:
                return
        else:  # REVIEWED - sanity test, with params re-write this logic/check can move there
            raise ValueError("Agents must be either IDU, NIDU, or ND")

    def _pca_interaction(self, relationship: Relationship, time):
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

        def influence(agent, partner):
            agent_opinion = agent.opinion
            partner_opinion = partner.opinion
            agent_influence = nx.closeness_centrality(self.get_Graph(), agent)
            partner_influence = nx.closeness_centrality(self.get_Graph(), partner)
            if agent_influence > partner_influence:
                partner.opinion = np.mean([agent.opinion, partner.opinion])
            elif agent_influence == partner_influence:
                pass
            else:
                agent.opinion = np.mean([agent.opinion, partner.opinion])
            if self.runRandom.random() < params.PCA_PrEP:
                if agent_opinion < params.opinion_threshold < agent.opinion:
                    self._initiate_PrEP(agent, time, force=True)
                elif partner_opinion < params.opinion_threshold < partner.opinion:
                    self._initiate_PrEP(partner, time, force=True)

        def knowledge_dissemination(partner):
            partner.awareness = True
            self.aware_agentSet.add_agent(partner)
            if (
                partner.opinion > params.opinion_threshold
                and self.runRandom.random() < params.PCA_PrEP
            ):
                self._initiate_PrEP(partner, time, force=True)

        def transmission_probability():
            if relationship._ID1.awareness and relationship._ID2.awareness:
                p = params.opinionTransmission
            else:
                p = params.knowledgeTransmission

            if num_acts == 1:
                p_total_transmission = p
            elif num_acts >= 1:
                p_total_transmission = 1.0 - binom.pmf(0, num_acts, p)
            else:
                p_total_transmission = 0
            return p_total_transmission

        acts_prob = self.runRandom.random()
        actsBin = 0
        current_p_value = 0.0

        while acts_prob > current_p_value:
            actsBin += 1
            current_p_value += params.interactionProb[relationship._rel_type][actsBin][
                "pvalue"
            ]

        minimum = params.interactionProb[relationship._rel_type][actsBin]["min"]
        maximum = params.interactionProb[relationship._rel_type][actsBin]["max"]
        if minimum == maximum:
            num_acts = minimum
        else:
            num_acts = self.runRandom.randrange(int(minimum), int(maximum))
        if num_acts < 1:
            return

        if relationship._ID1.awareness and not relationship._ID2.awareness:
            if self.runRandom.random() < transmission_probability():
                partner = relationship._ID2
                knowledge_dissemination(partner)
        elif not relationship._ID1.awareness and relationship._ID2.awareness:
            if self.runRandom.random() < transmission_probability():
                partner = relationship._ID2
                knowledge_dissemination(partner)
        elif relationship._ID1.awareness and relationship._ID2.awareness:
            if self.runRandom.random() < transmission_probability():
                influence(relationship._ID1, relationship._ID2)

    def _needle_transmission(self, agent: Agent, partner: Agent, time: int):
        """
        :Purpose:
            Simulate random transmission of HIV between two IDU agents
            through needle.\n
            Agent must by HIV+ and partner not.

        :Input:
            agents : int
            partner : int
            time : int
        :Output: -
        """

        assert agent._HIV_bool
        assert not partner._HIV_bool
        assert agent._DU == "IDU"
        assert partner._DU == "IDU"

        agent_race = agent._race
        agent_sex_type = agent._SO

        MEAN_N_ACTS = (
            params.DemographicParams[agent_race][agent_sex_type]["NUMSexActs"]
            * params.cal_NeedleActScaling
        )
        share_acts = poisson.rvs(MEAN_N_ACTS, size=1)[0]

        if agent._SNE_bool:  # Do they share a needle?
            p_UnsafeNeedleShare = 0.02  # no needle sharing
        else:  # they do share a needle
            # HIV+ ?
            if share_acts < 1:
                share_acts = 1

            p_UnsafeNeedleShare = (
                params.DemographicParams[agent_race][agent_sex_type]["NEEDLESH"]
                * params.safeNeedleExchangePrev
            )

        for n in range(share_acts):
            if self.runRandom.random() > p_UnsafeNeedleShare:
                share_acts -= 1

        if share_acts >= 1.0:
            p = agent.get_transmission_probability("NEEDLE")

            p_total_transmission: float
            if share_acts == 1:
                p_total_transmission = p
            else:
                p_total_transmission = 1.0 - binom.pmf(0, share_acts, p)

            if self.runRandom.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                self._become_HIV(partner, time)

    def _sex_transmission(self, time: int, rel: Relationship):
        """
        :Purpose:
            Simulate random transmission of HIV between two agents through Sex.
            Needed for all users. Sex is not possible in case the agent and
            assigned partner have incompatible Sex behavior. Given other logic,
            only one member of the relationship (the agent) has HIV at this time.

        :Input:
            time : int
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
        elif agent._HIV_bool or partner._HIV_bool:  # REVIEW: can't this just be else?
            # Define probabilities for unsafe sex
            # unprotected sex probabilities for primary partnerships
            MSexActs = (
                agent.get_number_of_sexActs(self.runRandom)
                * params.cal_SexualActScaling
            )
            T_sex_acts = int(poisson.rvs(MSexActs, size=1))

            num_int = rel._total_sex_acts
            # Get condom usage
            if params.condomUseType == "Race":
                p_UnsafeSafeSex1 = params.DemographicParams[agent._race][agent._SO][
                    "UNSAFESEX"
                ]
            else:
                p_UnsafeSafeSex1 = prob.unsafe_sex(num_int)

            # Reduction of risk acts between partners for condom usage
            # REVIEWED - what's with U_sex_acts1 and U_sex_acts2, U_sex_acts2 never seems to update - max to review
            U_sex_acts = T_sex_acts
            for n in range(U_sex_acts):
                if self.runRandom.random() < p_UnsafeSafeSex1:
                    U_sex_acts -= 1

            if U_sex_acts >= 1:
                # if agent HIV+
                rel._total_sex_acts += U_sex_acts
                ppAct = agent.get_transmission_probability("SEX")

                # Reduction of transmissibility for acts between partners for PrEP adherence
                if agent._PrEP_bool or partner._PrEP_bool:
                    if agent._PrEPresistance or partner._PrEPresistance:
                        pass

                    else:
                        if "Oral" in params.PrEP_type:  # params.PrEP_type == "Oral":
                            if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                                ppAct = ppAct * (1.0 - params.PrEP_AdhEffic)  # 0.04
                            else:
                                ppAct = ppAct * (1.0 - params.PrEP_NonAdhEffic)  # 0.24

                        elif "Inj" in params.PrEP_type:
                            ppActReduction = (
                                -1.0 * np.exp(-5.528636721 * partner._PrEP_load) + 1
                            )
                            if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                                ppAct = ppAct * (1.0 - ppActReduction)  # 0.04
                if partner.vaccine_bool or agent.vaccine_bool:
                    if params.vaccine_type == "HVTN702":
                        ppActPerc = 1 - np.exp(
                            -2.88
                            + 0.76 * (np.log10((partner.vaccine_time + 0.001) * 30))
                        )
                    elif params.vaccine_type == "RV144":
                        ppActPerc = 1 - np.exp(
                            -2.40 + 0.76 * (np.log10(partner.vaccine_time))
                        )
                    ppAct *= 1 - ppActPerc
                p_total_transmission: float
                if U_sex_acts == 1:
                    p_total_transmission = ppAct
                else:
                    p_total_transmission = 1.0 - binom.pmf(0, U_sex_acts, ppAct)

                if self.runRandom.random() < p_total_transmission:
                    # if agent HIV+ partner becomes HIV+
                    self._become_HIV(partner, time)

    def _become_HIV(self, agent: Agent, time: int):  # TODO rename
        """
        :Purpose:
            agent becomes HIV agent. Update all appropriate list and
            dictionaries.

        :Input:
            agent : int
            time : int
        """
        if not agent._HIV_bool:
            agent._HIV_bool = True
            agent._HIV_time = 1
            agent.vaccine_bool = False
            self.NewInfections.add_agent(agent)
            self.HIV_agentSet.add_agent(agent)
            if agent._PrEP_time > 0:
                if self.runRandom.random() < params.PrEP_resist:
                    agent._PrEPresistance = 1

        if agent._PrEP_bool:
            self._discont_PrEP(agent, time, force=True)

    def _enroll_treatment(self, time: int):
        """
        :Purpose:
            Enroll IDU agents in needle exchange

        :Input:
            time : int
        """
        print(("\n\n!!!!Engaginge treatment process: %d" % time))
        self.treatmentEnrolled = True
        for agent in self.All_agentSet.iter_agents():
            if self.runRandom.random() < params.treatmentCov and agent._DU == "IDU":
                agent._SNE_bool = True

    # REVIEWED this isn't used anywhere, but should be! _incarcerate makes things high risk and should reference this
    def _become_high_risk(self, agent: Agent, duration: int = None):

        if agent not in self.highrisk_agentsSet._members:
            self.highrisk_agentsSet.add_agent(agent)
        if not agent._everhighrisk_bool:
            self.NewHRrolls.add_agent(agent)

        agent._highrisk_bool = True
        agent._everhighrisk_bool = True

        if duration is not None:
            agent._highrisk_time = duration
        else:
            agent._highrisk_time = params.DemographicParams[agent._race][agent._SO][
                "HighRiskDuration"
            ]

    def _incarcerate(self, agent: Agent, time: int):
        """
        :Purpose:
            To incarcerate an agent or update their incarceration variables

        :Input:
            agent : int
            time : int

        """
        sex_type = agent._SO
        race_type = agent._race
        hiv_bool = agent._HIV_bool
        tested = agent._tested
        incar_t = agent._incar_time
        incar_bool = agent._incar_bool
        haart_bool = agent._HAART_bool
        recidivism = 1.0
        if agent._ever_incar_bool:
            recidivism = params.inc_Recidivism

        if incar_bool:
            agent._incar_time -= 1

            # get out if t=0
            if incar_t == 1:  # FREE AGENT
                self.incarcerated_agentSet.remove_agent(agent)
                self.NewIncarRelease.add_agent(agent)
                agent._incar_bool = False
                agent._ever_incar_bool = True
                if params.model == "Incar":
                    if (
                        not agent._highrisk_bool and params.flag_high_risk
                    ):  # If behavioral treatment on and agent HIV, ignore HR period.
                        if (
                            params.inc_treat_HRsex_beh
                            and hiv_bool
                            and (time >= params.inc_treatment_startdate)
                        ):
                            pass
                        else:  # Else, become high risk
                            self._become_high_risk(agent)

                            agent._mean_num_partners = (
                                agent._mean_num_partners + params.HR_partnerScale
                            )

                    if (
                        params.inc_treat_RIC
                        or params.inc_treat_HRsex_beh
                        or params.inc_treat_IDU_beh
                    ) and (time >= params.inc_treatment_startdate):
                        agent._incar_treatment_time = params.inc_treatment_dur

                    if hiv_bool:
                        if haart_bool:
                            if (
                                self.runRandom.random() > params.inc_ARTdisc
                            ):  # 12% remain surpressed
                                pass

                            else:
                                agent._HAART_bool = False
                                agent._HAART_adh = 0
                                self.Trt_ART_agentSet.remove_agent(agent)

                            # END FORCE

        elif self.runRandom.random() < params.DemographicParams[race_type][sex_type][
            "INCAR"
        ] * (1 + (hiv_bool * 4)) * params.cal_IncarP * (1 + recidivism):
            if agent._SO == "HF":
                jailDuration = (
                    prob.HF_jail_duration
                )  # TODO: make this into one dict for all SOs
            elif agent._SO == "HM":
                jailDuration = prob.HM_jail_duration

            durationBin = current_p_value = 0
            p = self.runRandom.random()
            while p > current_p_value:
                durationBin += 1
                current_p_value += jailDuration[durationBin]["p_value"]
            timestay = self.runRandom.randint(
                jailDuration[durationBin]["min"], jailDuration[durationBin]["max"]
            )

            if hiv_bool:
                if not tested:
                    if self.runRandom.random() < params.inc_PrisTestProb:
                        agent._tested = True
                else:  # Then tested and HIV, check to enroll in ART
                    if self.runRandom.random() < params.inc_ARTenroll:
                        tmp_rnd = self.runRandom.random()
                        HAART_ADH = params.inc_ARTadh
                        if tmp_rnd < HAART_ADH:
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
            for tmpA in agent._partners:
                if not tmpA._highrisk_bool:
                    if self.runRandom.random() < params.HR_proportion:
                        self._become_high_risk(tmpA)

                if params.flag_PrEP and (
                    params.PrEP_target_model == "incarcerated"
                    or params.PrEP_target_model == "incarcerated_high_risk"
                ):
                    # Atempt to put partner on prep if less than probability
                    if not tmpA._HIV_bool and not agent.vaccine_bool:
                        self._initiate_PrEP(tmpA, time)

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
                params.setting == "Scott"
            ):  # TODO fix this logic; should get partnerTraced and then lose it after
                # For each partner, determine if found by partner testing
                for ptnr in agent._partners:
                    if ptnr._HIV_bool and not ptnr._tested:
                        ptnr.partnerTraced = True
                        ptnr.traceTime = time + 1

        if not tested:
            test_prob = params.DemographicParams[race_type][sex_type]["HIVTEST"]

            # Rescale based on calibration param
            test_prob *= params.cal_TestFreq

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

        # Set HAART agents adherence at t=0 for all instanced HAART
        if time == 0 and agent_haart:
            if agent._HAART_adh == 0:
                HAART_ADH = params.DemographicParams[agent_race][agent_so]["HAARTadh"]
                if self.runRandom.random() < HAART_ADH:
                    adherence = 5
                else:
                    adherence = self.runRandom.randint(1, 4)

                # add to agent haart set
                agent._HAART_adh = adherence
                agent._HAART_time = time

        # Determine probability of HIV treatment
        if time >= 0 and agent._tested:
            # Go on HAART
            if not agent_haart and agent._HAART_time == 0:
                if (
                    self.runRandom.random()
                    < params.DemographicParams[agent_race][agent_so]["HAARTprev"]
                    * params.cal_ART_cov
                ):

                    HAART_ADH = params.DemographicParams[agent_race][agent_so][
                        "HAARTadh"
                    ]
                    if self.runRandom.random() < HAART_ADH:
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
                < params.DemographicParams[agent_race][agent_so]["HAARTdisc"]
            ):
                if not (agent._incar_treatment_time > 0 and params.inc_treat_RIC):
                    agent._HAART_bool = False
                    agent._HAART_adh = 0
                    agent._HAART_time = 0
                    self.Trt_ART_agentSet.remove_agent(agent)

    def _PrEP_eligible(
        self, agent: Agent, time: int
    ) -> bool:  # REVIEWED should this be in agent? self not used - move to agent
        eligible = False
        if len(params.PrEP_target_model & {"Allcomers", "Racial"}) > 0:
            eligible = True
        elif "CDCwomen" in params.PrEP_target_model:
            if agent._SO == "HF":
                for ptn in set(agent._relationships):
                    if ptn._ID1 == agent:
                        partner = ptn._ID2
                    else:
                        partner = ptn._ID1
                    if ptn._duration > 1:
                        if partner._DU == "IDU":
                            eligible = True
                            agent._PrEP_reason.append("IDU")
                        if partner._tested:
                            eligible = True
                            agent._PrEP_reason.append("HIV test")
                        if partner._MSMW:
                            eligible = True
                            agent._PrEP_reason.append("MSMW")
        elif "cdc_msm" in params.PrEP_target_model:
            if agent._SO == "MSM":
                for ptn in agent._relationships:
                    if ptn._ID1 == agent:
                        partner = ptn._ID2
                    else:
                        partner = ptn._ID1
                    if ptn._duration > 1:
                        if partner._tested or agent._mean_num_partners > 1:
                            eligible = True
                            break
        elif params.PrEP_target_model == "HighPN5":
            if agent._mean_num_partners >= 5:
                eligible = True
        elif (
            params.PrEP_target_model == "HighPN10"
        ):  # REVIEW is this used? if so, update to if x in y
            if agent._mean_num_partners >= 10:
                eligible = True
        elif params.PrEP_target_model == "MSM":
            if agent._SO == ("MSM" or "MTF"):
                eligible = True
        elif "RandomTrial" in params.PrEP_target_model and not params.flag_PCA:
            # If using random trial
            if time == 0:
                # if in init timestep 0, use agent set elligiblity
                eligible = agent.PrEP_eligible(
                    time
                )  # REVIEWED agent doesn't have this attribute, should this refer to the method? - use method (method needs to be moved, do that first)
            if time > 0:
                # else, false to not continue enrollment past random trial start
                eligible = False

        return eligible

    def _calc_PrEP_load(
        self, agent: Agent
    ):  # REVIEWED should this be in agent instead? self not used - move to agent
        """
        :Purpose:
            Determine load of PrEP concentration in agent.

        :Input:
            agent : agent()

        :Output:
            none
        """
        # N(t) = N0 (0.5)^(t/t_half)
        agent._PrEP_lastDose += 1
        if agent._PrEP_lastDose > 12:
            agent._PrEP_load = 0.0
        else:
            agent._PrEP_load = params.PrEP_peakLoad * (
                (0.5) ** (agent._PrEP_lastDose / (params.PrEP_halflife / 30))
            )

    def _discont_PrEP(self, agent: Agent, time: int, force: bool = False):
        # Agent must be on PrEP to discontinue PrEP
        assert agent._PrEP_bool

        # If force flag set, auto kick off prep.
        if force:
            self.Trt_PrEP_agentSet.remove_agent(agent)
            self.PrEPagents[agent._race][agent._SO] -= 1
            agent._PrEP_bool = False
            agent._PrEP_reason = []
            agent._PrEP_time = 0
        # else if agent is no longer enrolled on PrEP, increase time since last dose
        elif agent._PrEP_time > 0:
            agent._PrEP_time -= 1

        # else if agent is on PrEP, see if they should discontinue
        elif agent._PrEP_bool and agent._PrEP_time == 0:
            if (
                self.runRandom.random()
                < params.DemographicParams[agent._race][agent._SO]["PrEPdisc"]
            ):
                agent._PrEP_time = params.PrEP_falloutT
                self.Trt_PrEP_agentSet.remove_agent(agent)
                self.PrEPagents[agent._race][agent._SO] -= 1

                if "Oral" in params.PrEP_type:
                    agent._PrEP_bool = False
                    agent._PrEP_reason = []
            else:  # if not discontinue, see if its time for a new shot.
                if agent._PrEP_lastDose > 2:
                    agent._PrEP_lastDose = -1

        if "Inj" in params.PrEP_type:
            agent.update_PrEP_load()

    def vaccinate(self, ag: Agent, vax: str):
        ag.vaccine_bool = True
        ag.vaccine_type = vax
        ag.vaccine_time = 1

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
        initVaccine = params.init_with_vaccine
        if agent.vaccine_bool and not burn:
            agent.vaccine_time += 1
            if (
                params.flag_booster
                and agent.vaccine_time
                == params.DemographicParams[agent._race][agent._SO]["boosterInterval"]
                and self.runRandom.random()
                < params.DemographicParams[agent._race][agent._SO]["boosterProb"]
            ):
                self.vaccinate(agent, vaxType)

        elif time == params.vaccine_start:
            if initVaccine and burn:
                if (
                    self.runRandom.random()
                    < params.DemographicParams[agent._race][agent._SO]["vaccinePrev"]
                ):
                    self.vaccinate(agent, vaxType)
            elif not initVaccine and not burn:
                if (
                    self.runRandom.random()
                    < params.DemographicParams[agent._race][agent._SO]["vaccinePrev"]
                ):
                    self.vaccinate(agent, vaxType)

    def _initiate_PrEP(self, agent: Agent, time: int, force: bool = False):
        """
        :Purpose:
            Place agents onto PrEP treatment.
            PrEP treatment assumes that the agent knows his HIV+ status is negative.

        :Input:
            agent : Agent
            time : int
            force : default is `False`

        :Output:
            none
        """

        # REVIEWED _PrEP_time is initialized to zero, but then decremented to remove from PrEP - Sarah to review with Max
        def _enrollPrEP(self, agent: Agent):
            agent._PrEP_bool = True
            agent._PrEP_time = 0
            self.Trt_PrEP_agentSet.add_agent(agent)
            self.newPrEPagents.add_agent(agent)

            self.PrEPagents[agent._race][agent._SO] += 1
            if params.PrEP_target_model == "CDCwomen":
                if "IDU" in agent._PrEP_reason:
                    self.IDUprep += 1
                if "HIV test" in agent._PrEP_reason:
                    self.HIVprep += 1
                if "MSMW" in agent._PrEP_reason:
                    self.MSMWprep += 1

            tmp_rnd = self.runRandom.random()
            if params.setting == "AtlantaMSM":
                if (
                    tmp_rnd
                    < params.DemographicParams[agent._race][agent._SO]["PrEPadh"]
                ):
                    agent._PrEP_adh = 1
                else:
                    agent._PrEP_adh = 0
            else:
                if tmp_rnd < params.PrEP_Adherence:
                    agent._PrEP_adh = 1
                else:
                    agent._PrEP_adh = 0

            # set PrEP load and dosestep for PCK
            if "Inj" in params.PrEP_type:
                agent._PrEP_load = params.PrEP_peakLoad
                agent._PrEP_lastDose = 0
                # TODO: make this work for vaccines, all prep types
                if self.runRandom.random() < params.LAI_chance:
                    self.LAI_agentSet.add_agent(agent)
                    agent.PrEP_type = "LAI"
                else:
                    self.oralPrEP_agentSet.add_agent(agent)
                    agent.PrEP_type = "Oral"

        # agent must exist
        assert agent is not None

        # Check valid input

        # Prep only valid for agents not on prep and are HIV negative
        if agent._PrEP_bool or agent._HIV_bool:
            return

        # Determine probability of HIV treatment
        agent_race = agent._race

        if force:
            _enrollPrEP(self, agent)
        else:
            if "Racial" in params.PrEP_target_model:
                numPrEP_agents = sum(self.PrEPagents[agent_race].values())
            else:
                numPrEP_agents = self.Trt_PrEP_agentSet.num_members()

            if params.PrEP_target_model == "Clinical":
                target_PrEP_population = (
                    self.All_agentSet.num_members() - self.HIV_agentSet.num_members()
                )
                target_PrEP = target_PrEP_population * params.PrEP_Target

            elif (
                params.PrEP_target_model == "Incar"
                or params.PrEP_target_model == "IncarHR"
            ):
                if self.runRandom.random() < params.PrEP_Target:
                    _enrollPrEP(self, agent)
                return None
            elif "Racial" in params.PrEP_target_model:
                all_HIV_agents = set(self.All_agentSet._subset["HIV"]._members)

                all_race = set(
                    self.All_agentSet._subset["Race"]._subset[agent._race]._members
                )
                HIV_agents = len(all_HIV_agents & all_race)
                target_PrEP = (
                    int(
                        self.All_agentSet._subset["Race"]
                        ._subset[agent._race]
                        .num_members()
                    )
                    - HIV_agents
                ) * params.DemographicParams[agent._race][agent._SO]["PrEP_coverage"]

            else:
                target_PrEP = int(
                    (
                        self.All_agentSet.num_members()
                        - self.All_agentSet._subset["HIV"].num_members()
                    )
                    * params.PrEP_Target
                )

            if params.PrEP_clinic_cat == "Racial" and agent_race == "BLACK":
                if self.runRandom.random() < params.PrEP_Target:
                    _enrollPrEP(self, agent)
            elif (
                params.PrEP_target_model == "Incar"
                or params.PrEP_target_model == "IncarHR"
            ):
                if self.runRandom.random() < params.PrEP_Target:
                    _enrollPrEP(self, agent)
            elif (
                numPrEP_agents < target_PrEP
                and time >= params.PrEP_startT
                and agent.PrEP_eligible(time)
            ):
                _enrollPrEP(self, agent)

    def _get_clinic_agent(
        self, clinicBin: str, eligiblePool: Sequence[Agent]
    ) -> Optional[Agent]:
        pMatch = 0.0
        RN = self.runRandom.random()
        for i in range(6):  # 6 is exclusive
            if RN < pMatch:
                break
            else:
                pMatch += params.clinicAgents[clinicBin][i]["Prob"]
                minNum = params.clinicAgents[clinicBin][i]["min"]
                maxNum = params.clinicAgents[clinicBin][i]["max"]

        # try twice to match
        for i in range(1, 3):
            randomK_sample = self.runRandom.sample(
                eligiblePool, params.cal_ptnrSampleDepth
            )
            eligibleK_Pool = [
                ag
                for ag in randomK_sample
                if (
                    (ag._mean_num_partners >= minNum)
                    and (ag._mean_num_partners <= maxNum)
                )
            ]
            if eligibleK_Pool:
                return self.runRandom.choice(eligibleK_Pool)
            else:
                print(
                    (
                        "Looking for agent with min:%d and max %d failed %d times"
                        % (minNum, maxNum, i)
                    )
                )

        print("No suitable PrEP agent")
        return None

    def _progress_to_AIDS(self, agent: Agent):
        """
        :Purpose:
            Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        if not agent._HIV_bool:
            raise ValueError("AIDS only valid for HIV agents!agent:%s" % str(agent._ID))

        # if agent not in self.AIDS_agents:
        if not agent._HAART_bool:
            adherenceStat = agent._HAART_adh
            p = prob.adherence_prob(adherenceStat)

            if self.runRandom.random() < p * params.cal_ProgAIDS:
                agent._AIDS_bool = True
                self.HIV_AIDS_agentSet.add_agent(agent)

    def _die_and_replace(self, time: int, reported: bool = True):

        """
        :Purpose:
            Let agents die and replace the dead agent with a new agent randomly.
        """
        for agent in self.All_agentSet._members:

            # agent incarcerated, don't evaluate for death
            if agent._incar_bool:
                pass

            # death rate per 1 person-month
            p = prob.get_death_rate(
                agent._HIV_bool, agent._AIDS_bool, agent._race, agent._HAART_adh
            )

            if self.runRandom.random() < p:
                self.deathSet.append(agent)

                # End all existing relationships
                for rel in agent._relationships:
                    rel.progress(forceKill=True)

                    self.Relationships.remove(rel)

                # Remove agent node and edges from network graph
                self.get_Graph().remove_node(agent)

                # Remove agent from agent class and sub-sets
                self.All_agentSet.remove_agent(agent)

                agent_cl = self._return_new_Agent_class(agent._race, agent._SO)
                self.create_agent(agent_cl, agent._race)
                self.get_Graph().add_node(agent_cl)
