#!/usr/bin/env python3
# encoding: utf-8

# Imports
import random
from random import Random
from typing import Dict, Any, List, Sequence, Optional
import os

import numpy as np  # type: ignore
from scipy.stats import binom  # type: ignore
from scipy.stats import poisson  # type: ignore
import networkx as nx  # type: ignore

try:
    from .agent import Agent_set, Agent, Relationship
except ImportError as e:
    raise ImportError("Can't import network_graph_tools! %s" % str(e))

try:
    from .network_graph_tools import NetworkClass
except ImportError as e:
    raise ImportError("Can't import network_graph_tools! %s" % str(e))

try:
    from .analysis_output import initiate_ResultDict, print_stats
except ImportError as e:
    raise ImportError("Can't import analysis_output! %s" % str(e))

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

        :py:class:`NetworkClass` : Inherited
        :py:class:`PopulationClass` : Inherited

    :Attributes:
        :py:attr:`tmax` : int
            Number of time steps simulated.

        :py:attr:`CleanSyringeUsers` : list

        :py:attr:`SEPAgents` : dict
            Dictionary of users who participated in a
            syringe exchange program (SEP) {agent:time}.

        :py:attr:`DrugTreatmentAgents` : dict
            Dictionary of users who underwent drug
            treatment {agent:time}.

        :py:attr:`TestedAgents` : list
            List of agents who get tested for HIV every time step.

        :py:attr:`tmp_Agents` : dict
            Changes resulting from parsing through the agents
            and applying the update rules are stored
            in :py:attr:`tmp_agent_dict`.

        All attributes from :py:class:`NetworkClass` \n

        All attributes from :py:class:`PopulationClass`

    :Methods:
        :py:meth:`_update_population` \n
        :py:meth:`_needle_transmission` \n
        :py:meth:`_sex_transmission` \n
        :py:meth:`_update_IDU` \n
        :py:meth:`_update_AllAgents` \n
        :py:meth:`_initiate_HAART` \n
        :py:meth:`_update_AllAgents` \n
        :py:meth:`run` \n
        All methods from :py:class:`NetworkClass` \n
        All methods from :py:class:`PopulationClass`
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
        parameter_dict: Dict[str, Any],
        runseed: int,
        popseed: int,
        netseed: int,
        network_type: str = "",
    ):
        """ Initialize HIVModel object """
        # Ensure param variable is are defined. For backwards compatibility with params.py files
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

        self.current_dir = os.getcwd()
        print("=== Begin Initialization Protocol ===\n")

        NetworkClass.__init__(
            self,
            N=N,
            network_type=network_type,
            popSeed=self.popseed,
            netSeed=self.netseed,
        )

        # keep track of current time step globally for dynnetwork report
        self.TimeStep = 0
        self.totalIncarcerated = 0

        print("\n\tCreating lists")
        # Other lists / dictionaries

        self.NewInfections = Agent_set("NewInfections")
        self.NewDiagnosis = Agent_set("NewDiagnosis")
        self.NewIncarRelease = Agent_set("NewIncarRelease")
        self.NewHRrolls = Agent_set("NewHRrolls")

        self.Acute_agents: List[Agent] = []
        self.Transmit_from_agents: List[Agent] = []
        self.Transmit_to_agents: List[Agent] = []

        self.totalDiagnosis = 0
        self.treatmentEnrolled = False
        self.num_Deaths = self._reset_death_count()

        self.ResultDict = initiate_ResultDict()
        self.newPrEPagents = Agent_set("NewPrEPagents")
        self.newPrEPenrolls = 0
        self.IDUprep = 0
        self.HIVprep = 0
        self.MSMWprep = 0
        # Set seed format. 0: pure random, -1: Stepwise from 1 to nRuns, else: fixed value

        print(("\tRun seed was set to:", runseed))
        self.runRandom = Random(
            runseed
        )  #  - what if self.runseed != runseed (if 0 passed)
        random.seed(self.runseed)
        np.random.seed(self.runseed)
        print(("\tFIRST RANDOM CALL %d" % random.randint(0, 100)))

        print("\tResetting death count")
        self._reset_death_count()  # Number of death

        print("\tCreating network graph")
        self.create_graph_from_agents(
            self.All_agentSet
        )  # REVIEWED redundant with NetworkClass init? - review with max, logic feels scattered as NetworkClass also intializes a graph

        print("\n === Initialization Protocol Finished ===")

    def run(self, dir_prefix="Results"):
        """
        Core of the model:
            1. Prints networkReport for first agents.
            2. Makes agents become HIV (used for current key_time tracking for acute)
            3. Loops over all time steps
                a. _update AllAgents()
                b. _reset_death_counts()
                c. _ self._die_and_replace()
                d. self._update_population()
                e. self._reset_partner_count()
        """

        def getStats(t):
            print_stats(
                self,
                self.runseed,
                t,
                self.All_agentSet,
                self.HIV_agentSet,
                self.incarcerated_agentSet,
                self.Trt_PrEP_agentSet,
                self.NewInfections,
                self.NewDiagnosis,
                self.num_Deaths,
                self.ResultDict,
                self.Relationships,
                self.NewHRrolls,
                self.NewIncarRelease,
                self.deathSet,
            )

        def print_components(t):
            name = "componentReport_ALL"
            compReport = open("results/" + name + ".txt", "a")
            components = sorted(
                nx.connected_component_subgraphs(self.get_Graph()),
                key=len,
                reverse=True,
            )
            compID = 0
            for comp in components:
                totN = nhiv = ntrtmt = ntrthiv = nprep = PrEP_ever_HIV = 0
                for ag in comp.nodes():
                    totN += 1
                    if ag._HIV_bool:
                        nhiv += 1
                        if ag._treatment_bool:
                            ntrthiv += 1
                        if ag._PrEP_ever_bool:
                            PrEP_ever_HIV += 1
                    elif ag._treatment_bool:
                        ntrtmt += 1
                        if ag._PrEP_bool:
                            nprep += 1
                compReport.write(
                    "{rseed}\t{pseed}\t{nseed}\t{t}\t{compID}\t{totalN}\t{Nhiv}\t{Ntrtmt}\t{Nprep}\t{NtrtHIV}\t{NprepHIV}\n".format(
                        rseed=self.runseed,
                        pseed=self.popseed,
                        nseed=self.netseed,
                        t=t,
                        compID=compID,
                        totalN=totN,
                        Nhiv=nhiv,
                        Ntrtmt=ntrtmt,
                        Nprep=nprep,
                        NtrtHIV=ntrthiv,
                        NprepHIV=PrEP_ever_HIV,
                    )
                )

                compID += 1
            compReport.close()

        def burnSimulation(burnDuration: int):
            print(
                ("\n === Burn Initiated for {} timesteps ===".format(burnDuration + 1))
            )
            for t in range(0, burnDuration + 1):
                self._update_AllAgents(t, burn=True)

                if params.flag_DandR:
                    self._die_and_replace(t)
            print(("\tBurn Cuml Inc:\t{}".format(self.NewInfections.num_members())))
            self.NewInfections.clear_set()
            self.NewDiagnosis.clear_set()
            self.NewHRrolls.clear_set()
            self.NewIncarRelease.clear_set()
            self.newPrEPagents.clear_set()
            self.newPrEPenrolls = 0
            self.IDUprep = 0
            self.HIVprep = 0
            self.MSMWprep = 0

            self._reset_death_count()
            print(" === Simulation Burn Complete ===")

        burnSimulation(params.burnDuration)

        print("\n === Begin Simulation Run ===")
        if params.drawFigures:
            nNodes = self.G.number_of_nodes()
            self.visualize_network(
                coloring=params.drawFigureColor,
                node_size=5000.0 / nNodes,
                curtime=0,
                iterations=10,
                label="Seed" + str(self.runseed),
            )
        if params.calcComponentStats:
            print_components(0)

        self.cumInfT = 0
        self.cumInfW = 0
        self.cumInfB = 0

        def makeAgentZero(numPartners: int):
            firstHIV = self.runRandom.choice(self.DU_IDU_agentSet._members)
            i = 0
            while i <= numPartners:
                self.update_partner_assignments(
                    10000.0, self.get_Graph(), agent=firstHIV
                )
                i += 1
            self._become_HIV(firstHIV, 0)

        print("\t===! Start Main Loop !===")

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
                (
                    "\tSTARTING HIV count:{}\tTotal Incarcerated:{}\tHR+:{}\tPrEP:{}".format(
                        self.HIV_agentSet.num_members(),
                        self.incarcerated_agentSet.num_members(),
                        self.highrisk_agentsSet.num_members(),
                        self.Trt_PrEP_agentSet.num_members(),
                    )
                )
            )
            print(
                (
                    "Trt:{trt} \t OAT:{oat} \t NAL:{nal}".format(
                        trt=self.treatment_agentSet.num_members(),
                        oat=len(
                            [
                                a
                                for a in self.treatment_agentSet._members
                                if a._OAT_bool is True
                            ]
                        ),
                        nal=len(
                            [
                                a
                                for a in self.treatment_agentSet._members
                                if a._naltrex_bool is True
                            ]
                        ),
                    )
                )
            )
            self.TimeStep = t

            self._update_AllAgents(t)

            getStats(t)

            self._reset_death_count()

            if params.flag_DandR:
                self._die_and_replace(t)

            print(("Number of relationships: %d" % self.Relationships.num_members()))
            self.All_agentSet.print_subsets()

            newInfB = len(
                [tmpA for tmpA in self.NewInfections._members if tmpA._race == "BLACK"]
            )
            newInfW = len(
                [tmpA for tmpA in self.NewInfections._members if tmpA._race == "WHITE"]
            )
            newInfT = len(self.NewInfections._members)
            self.cumInfB += newInfB
            self.cumInfW += newInfW
            self.cumInfT += newInfT

            self.totalDiagnosis += len(self.NewDiagnosis._members)
            if (
                self.totalDiagnosis > params.initTreatment
                and not self.treatmentEnrolled
            ):
                self._enroll_treatment(t)

            self.NewInfections.clear_set()
            self.NewDiagnosis.clear_set()
            self.NewHRrolls.clear_set()
            self.NewIncarRelease.clear_set()
            # self.num_Deaths = {} # REVIEWED - not expected to be an empty dict (expects sub dicts) - should this be _reset_death_count? check if something else breaks (num deaths)
            prepReport = open("results/PrEPReport.txt", "a")
            prepReport.write(
                f"{self.runseed}\t{self.TimeStep}\t{self.newPrEPenrolls}\t{self.IDUprep}\t{self.HIVprep}\t{self.MSMWprep}\n"
            )
            prepReport.close()
            self.newPrEPagents.clear_set()
            self.newPrEPenrolls = 0
            self.IDUprep = 0
            self.HIVprep = 0
            self.MSMWprep = 0
            # If set to draw the edge list, print list at each timestep
            if params.drawEdgeList and t % params.intermPrintFreq == 0:
                print("Drawing network edge list to file")
                fh = open("results/network/Edgelist_t{}.txt".format(t), "wb")
                self.write_G_edgelist(fh)
                fh.close()

            print((t % params.intermPrintFreq))
            if t % params.intermPrintFreq == 0:
                if params.calcNetworkStats:
                    self.write_network_stats(t=t)
                if params.calcComponentStats:
                    print_components(t)

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

        self.Acute_agents = []
        self.Transmit_from_agents = []
        self.Transmit_to_agents = []

        for rel in self.Relationships._members:
            # If in burn, ignore interactions
            if burn:
                pass
            else:
                self._agents_interact(rel._ID1, rel._ID2, time, rel)

            # If static network, ignore relationship progression
            if params.flag_staticN:
                pass
            else:
                if rel.progress():
                    g = self.get_Graph()
                    if g.has_edge(rel._ID1, rel._ID2):
                        g.remove_edge(rel._ID1, rel._ID2)

                    self.Relationships.remove_agent(rel)
                    del rel

        if params.flag_HR:
            for tmpA in self.highrisk_agentsSet.iter_agents():
                if tmpA._highrisk_time > 0:
                    tmpA._highrisk_time -= 1
                    if (
                        tmpA._SO == "HM"
                        and params.flag_PrEP
                        and (
                            params.PrEP_target_model == "HR"
                            or params.PrEP_target_model == "IncarHR"
                        )
                    ):
                        for part in tmpA._partners:
                            if part._HIV_bool or part.vaccine_time >= 1:
                                pass
                            else:
                                self._initiate_PrEP(part, time)
                else:
                    self.highrisk_agentsSet.remove_agent(tmpA)
                    tmpA._highrisk_bool = False

                    if params.model == "Incar":
                        if tmpA._SO == "HM":
                            tmpA._mean_num_partners -= params.HR_partnerScale
                        elif tmpA._SO == "HF":
                            tmpA._mean_num_partners -= params.HR_partnerScale

        for agent in self.All_agentSet.iter_agents():

            agent_drug_type = agent._DU
            agent_HIV_status = agent._HIV_bool

            agent._timeAlive += 1

            if params.flag_incar:  # and not burn:
                self._incarcerate(agent, time)

            if agent._MSMW and self.runRandom.random() < params.HIV_MSMW:
                self._become_HIV(agent, 0)

            if agent_HIV_status:
                # If in burnin, ignore HIV
                if burn:
                    if agent._incar_treatment_time >= 1:
                        agent._incar_treatment_time -= 1

                else:
                    self._HIVtest(agent, time)
                    self._progress_to_AIDS(agent, agent_drug_type)

                    if params.flag_ART:
                        self._initiate_HAART(agent, time)
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
                        elif (
                            self._PrEP_eligible(agent, time)
                            and not agent._PrEP_bool
                            and agent.vaccine_time == 0
                        ):
                            self._initiate_PrEP(agent, time)
                if "Vaccine" in params.PrEP_type and not agent._PrEP_bool and not burn:
                    self.initiate_Vaccine(agent, time, vaxType=params.vaccine_type)

        if params.flag_PrEP and time >= params.PrEP_startT:
            if params.PrEP_target_model == "Clinical":
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
                params.PrEP_target_model == "RandomTrial" and time == params.PrEP_startT
            ):
                print("Starting random trial")
                components = sorted(
                    nx.connected_component_subgraphs(self.G), key=len, reverse=True
                )
                totNods = 0
                for comp in components:
                    totNods += comp.number_of_nodes()
                    if self.runRandom.random() < 0.5:
                        # Component selected as treatment pod!
                        for ag in comp.nodes():
                            if (ag._HIV_bool is False) and (ag._PrEP_bool is False):
                                ag._treatment_bool = True
                                if (
                                    self.runRandom.random() < params.PrEP_Target
                                    and agent.vaccine_time == 0
                                ):
                                    self._initiate_PrEP(ag, time, force=True)
                print(("Total agents in trial: ", totNods))

    def _agents_interact(
        self, agent_1: Agent, agent_2: Agent, time: int, rel: Relationship
    ):  # REVIEWED if rel includes agents, why pass the agents at all?, also only thing from self that is used is the random number generator, maybe pass this and move this to ABM_partnering - remove agents, use rel to get agents and update the logic below
        """
        :Purpose:
            Let IDU agent interact with a partner.
            Update IDU agents:
                1 - determine transition type
                2 - Injection rules
                3 - Sex rules
                4 - HIV transmission

        :Input:
            agent : int

            partner : int

            time : int

        Output:
            none

        """
        # print agent
        partner_HIV_status = agent_2._HIV_bool
        agent_HIV_status = agent_1._HIV_bool
        agent_incar = agent_1._incar_bool
        partner_incar = agent_2._incar_bool

        eligible = False

        # If either agent is incarcerated, skip their interaction
        if agent_incar or partner_incar:
            return
        # Else if neither agent is HIV (shouldn't be possible), skip their interaction to save computation time
        elif not agent_HIV_status and not partner_HIV_status:
            return
        elif agent_HIV_status:  # If agent_1 is HIV
            if (
                partner_HIV_status
            ):  # If agent_1 and agent_2 are both HIV, skip interaction
                return
            else:  # Agent is HIV, partner is succept
                agent = agent_1
                partner = agent_2
                eligible = True
        elif (
            partner_HIV_status
        ):  # If agent_2 is HIV and we have tested both HIV +/-, agent_2 is HIV, agent_1 is succept
            agent = agent_2
            partner = agent_1
            eligible = True

        if eligible:
            partner_drug_type = partner._DU
            agent_drug_type = agent._DU
            partner_sex_type = partner._SO
            agent_sex_type = agent._SO
            partner_HIV_status = partner._HIV_bool
            agent_HIV_status = agent._HIV_bool
            agent_incar = agent._incar_bool
            partner_incar = partner._incar_bool
            if partner_drug_type == "IDU" and agent_drug_type == "IDU":
                # Injection is possible
                # If agent is on post incar HR treatment to prevent IDU behavior, pass IUD infections
                if agent._incar_treatment_time > 0 and params.inc_treat_IDU_beh:
                    pass

                elif sex_possible(agent_sex_type, partner_sex_type):
                    # Sex is possible
                    rv = self.runRandom.random()
                    if rv < 0.25:  # Needle only (60%)
                        self._needle_transmission(agent, partner, time)
                    else:  # Both sex and needle (20%)
                        self._needle_transmission(agent, partner, time)
                        self._sex_transmission(
                            agent, partner, time, rel
                        )  # , num_interactions)
                else:
                    # Sex not possible, needle only
                    self._needle_transmission(agent, partner, time)

            elif partner_drug_type in ["NIDU", "NDU"] or agent_drug_type in [
                "NIDU",
                "NDU",
            ]:
                if sex_possible(agent_sex_type, partner_sex_type):
                    self._sex_transmission(agent, partner, time, rel)
                else:
                    return
            else:
                raise ValueError("Agents must be either IDU, NIDU, or ND")

    def get_acute_status(
        self, agent: Agent, time: int
    ) -> bool:  # REVIWED should this be in agent? self not used - move this
        """
        :Purpose:
            Simulate random transmission of HIV between two IDU agents
            through needle.\n
            Needed in _update_IDUand
        :Input:
            agents : int
            partner : int
        time : int
        :Output: -
        """
        acuteTimePeriod = 2
        hiv_t = agent._HIV_time

        if hiv_t <= acuteTimePeriod and hiv_t > 0:
            return True
        else:
            return False

    def get_transmission_probability(
        self, agent: Agent, interaction: str
    ) -> float:  # REVIEWED should this be in agent, self not (really) used - move to agent
        """ Decriptor
            :Purpose:
            Determines the probability of a transmission event based on type. Determines if act is needle/sexual,

            :Input:
                N : int
                Number of agents. Default: 1000
                tmax: int
                Number of simulation steps (years).

                :py:class:`NetworkClass` : Inherited
                :py:class:`PopulationClass` : Inherited

            :Attributes:
                :py:attr:`tmax` : int
                    Number of time steps simulated.
                """

        sex_type = agent._SO
        race_type = agent._race
        tested = agent._tested
        onHAART = agent._HAART_bool
        agentAdherence = agent._HAART_adh

        # Logic for if needle or sex type interaction
        p: float
        if interaction == "NEEDLE":
            p = params.TransmissionProbabilities["NEEDLE"][str(agentAdherence)]
        elif interaction == "SEX":
            p = params.TransmissionProbabilities["SEX"][sex_type][str(agentAdherence)]

        isAcute = self.get_acute_status(agent, 0)

        # Scaling parameter for acute HIV infections
        if isAcute:
            p = p * params.cal_AcuteScaling

        # Scaling parameter for positively identified HIV agents
        if tested:
            p = p * (1 - params.cal_RR_Dx)

        # Tuning parameter for ART efficiency
        if onHAART:  # self.AdherenceAgents[agent] > 0:
            p = p * params.cal_RR_HAART

        # Racial calibration parameter to attain proper race incidence disparity
        if race_type == "BLACK":
            p = p * params.cal_raceXmission

        # Scaling parameter for per act transmission.
        p = p * params.cal_pXmissionScaling

        return p

    def _needle_transmission(self, agent: Agent, partner: Agent, time: int):
        """
        :Purpose:
            Simulate random transmission of HIV between two IDU agents
            through needle.\n
            Needed in _update_IDUand
        :Input:
            agents : int
            partner : int
        time : int
        :Output: -
        """
        # Param to scale number of partners

        # both must be IDU
        partner_drug_type = partner._DU
        agent_drug_type = agent._DU
        agent_race = agent._race
        agent_sex_type = agent._SO
        Race_Agent = agent._race
        Type_agent = agent._SO

        if not (partner_drug_type == "IDU" and agent_drug_type == "IDU"):
            raise ValueError(
                "To share a needle both agents must be IDU!%s %s"
                % (str(agent_drug_type), str(partner_drug_type))
            )

        # Do they share a needle?
        SEPstat = agent._SNE_bool

        isAcute = self.get_acute_status(agent, time)
        HIV_agent = agent._HIV_bool
        HIV_partner = partner._HIV_bool
        MEAN_N_ACTS = (
            params.DemographicParams[Race_Agent][Type_agent]["NUMSexActs"]
            * params.cal_NeedleActScaling
        )
        share_acts = poisson.rvs(MEAN_N_ACTS, size=1)[0]

        if SEPstat:
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

        if HIV_agent == 1 and HIV_partner == 0 and share_acts >= 1.0:
            p = self.get_transmission_probability(agent, "NEEDLE")

            p_total_transmission: float
            if share_acts == 1:
                p_total_transmission = p
            else:
                p_total_transmission = 1.0 - binom.pmf(0, share_acts, p)

            if self.runRandom.random() < p_total_transmission:
                # if agent HIV+ partner becomes HIV+
                self._become_HIV(partner, time)
                self.Transmit_from_agents.append(agent)
                self.Transmit_to_agents.append(partner)
                if isAcute:
                    self.Acute_agents.append(agent)
                else:
                    pass

    def _sex_transmission(
        self, agent: Agent, partner: Agent, time: int, rel: Relationship
    ):  # REVIEWED rel contains agent and partner, why pass agent and partner? - just pass rel - update logic

        """
        :Purpose:
            Simulate random transmission of HIV between two agents through Sex.
            Needed for all users. Sex is not possible in case the agent and
            assigned partner have incompatible Sex behavior.

        :Input:
            agents : int
            partner : int
            time : int
        number_of_interaction : int

        :Output:
            none
        """
        # Double check: Sex possible?
        Type_agent = agent._SO
        Type_partner = partner._SO
        if not sex_possible(Type_agent, Type_partner):
            raise ValueError(
                "Sex must be possible! %s %s" % (str(Type_agent), str(Type_partner))
            )

        # HIV status of agent and partner
        # Everything from here is only run if one of them is HIV+

        HIVstatus_Agent = agent._HIV_bool
        HIVstatus_Partner = partner._HIV_bool
        Race_Agent = agent._race
        isAcute = self.get_acute_status(agent, time)

        if HIVstatus_Partner:
            pass
        if HIVstatus_Agent and HIVstatus_Partner:
            return
        elif HIVstatus_Agent or HIVstatus_Partner:
            # Define probabilities for unsafe sex
            # unprotected sex probabilities for primary partnerships
            MSexActs = self._get_number_of_sexActs(agent) * params.cal_SexualActScaling
            T_sex_acts1 = int(poisson.rvs(MSexActs, size=1))

            num_int = rel._total_sex_acts
            # Get condom usage
            if params.condomUseType == "Race":
                p_UnsafeSafeSex1 = params.DemographicParams[Race_Agent][Type_agent][
                    "UNSAFESEX"
                ]
            else:
                p_UnsafeSafeSex1 = prob.unsafe_sex(num_int)

            # Reduction of risk acts between partners for condom usage
            # REVIEWED - what's with U_sex_acts1 and U_sex_acts2, U_sex_acts2 never seems to update - max to review
            U_sex_acts1 = T_sex_acts1
            for n in range(U_sex_acts1):
                if self.runRandom.random() < p_UnsafeSafeSex1:
                    U_sex_acts1 -= 1

            U_sex_acts2 = U_sex_acts1

            if U_sex_acts2 >= 1:
                # if agent HIV+
                rel._total_sex_acts += U_sex_acts2
                if HIVstatus_Agent == 1 or HIVstatus_Partner == 1:
                    ppAct = self.get_transmission_probability(agent, "SEX")

                    # Reduction of transmissibility for acts between partners for PrEP adherence
                    if agent._PrEP_bool or partner._PrEP_bool:
                        if agent._PrEPresistance or partner._PrEPresistance:
                            pass

                        else:
                            if (
                                "Oral" in params.PrEP_type
                            ):  # params.PrEP_type == "Oral":
                                if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                                    ppAct = ppAct * (1.0 - params.PrEP_AdhEffic)  # 0.04
                                else:
                                    ppAct = ppAct * (
                                        1.0 - params.PrEP_NonAdhEffic
                                    )  # 0.24

                            elif "Inj" in params.PrEP_type:
                                ppActReduction = (
                                    -1.0 * np.exp(-5.528636721 * partner._PrEP_load) + 1
                                )
                                if agent._PrEP_adh == 1 or partner._PrEP_adh == 1:
                                    ppAct = ppAct * (1.0 - ppActReduction)  # 0.04
                    if partner.vaccine_time >= 1:
                        if params.vaccine_type == "HVTN702":
                            ppActReduction = 1 - np.exp(
                                -2.88
                                + 0.76 * (np.log(partner.vaccine_time + 0.001 * 30))
                            )

                        elif params.vaccine_type == "RV144":
                            ppActReduction = 1 - np.exp(-2.40 + 0.76 * (np.log(partner.vaccine_time)))
                        ppAct *= 1 - ppActReduction


                    p_total_transmission: float
                    if U_sex_acts2 == 1:
                        p_total_transmission = ppAct
                    else:
                        p_total_transmission = 1.0 - binom.pmf(0, U_sex_acts1, ppAct)

                    if self.runRandom.random() < p_total_transmission:

                        # if agent HIV+ partner becomes HIV+
                        self.Transmit_from_agents.append(agent)
                        self.Transmit_to_agents.append(partner)
                        self._become_HIV(partner, time)
                        if isAcute:
                            self.Acute_agents.append(agent)
                        else:
                            pass
            else:
                return

    def _get_number_of_sexActs(self, agent: Agent):
        """
        :Purpose:
            agent becomes HIV agent. Update all appropriate list and
            dictionaries.

        :Input:
            agent : int

        """
        # 1 time per year 96 1.9 29 0.9 67 3.4
        # 2–5 times per year 428 8.2 184 5.8 244 12.2
        # 6–11 times per year 328 6.3 183 5.7 145 7.3
        # 12–23 times per year 376 7.2 251 7.9 125 6.3
        # 24–35 times per year 1,551 29.9 648 20.3 903 45.3
        # 36–51 times per year 1,037 20.0 668 20.9 369 18.5
        # 52–155 times per year 644 12.4 605 18.9 39 2.0
        # >156 times per year 733 14.1 631 19.7 102 5.1
        rv = self.runRandom.random()
        pMatch = 0.0
        i = 0

        while True:
            i += 1
            pMatch += params.sexualFrequency[i]["p_value"]
            if rv <= pMatch:
                minSA = params.sexualFrequency[i]["min"]
                maxSA = params.sexualFrequency[i]["max"]
                return self.runRandom.randrange(minSA, maxSA, 1)
            if i == 5:
                break

    def _become_HIV(self, agent: Agent, time: int):
        """
        :Purpose:
            agent becomes HIV agent. Update all appropriate list and
            dictionaries.

        :Input:
            agent : int
            time : int

        """
        if agent._HIV_bool:
            pass
        else:
            agent._HIV_bool = True
            agent._HIV_time = 1
            agent.vaccine_time = 0
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
            Enroll agents with HIV in ART

        :Input:
            time : int

        """
        print(("\n\n!!!!Engaginge treatment process: %d" % time))
        self.treatmentEnrolled = True
        for agent in self.All_agentSet.iter_agents():
            if self.runRandom.random() < params.treatmentCov and agent._DU == "IDU":
                agent._SNE_bool = True

    def _becomeHighRisk(self, agent: Agent, HRtype: str = "", duration: int = None):

        if agent not in self.highrisk_agentsSet._members:
            self.highrisk_agentsSet.add_agent(agent)
        if not agent._everhighrisk_bool:
            self.NewHRrolls.add_agent(agent)

        agent._highrisk_bool = True
        agent._everhighrisk_bool = True
        agent._highrisk_type = HRtype

        if duration is not None:
            agent._highrisk_time = duration
        else:
            agent._highrisk_time = params.HR_M_dur

    def _incarcerate(self, agent: Agent, time: int):
        """
        :Purpose:
            To incarcerate agents

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
                        not agent._highrisk_bool and params.flag_HR
                    ):  # If behavioral treatment on and agent HIV, ignore HR period.
                        if (
                            params.inc_treat_HRsex_beh
                            and hiv_bool
                            and (time >= params.inc_treatment_startdate)
                        ):
                            pass
                        else:  # Else, become high risk
                            self.highrisk_agentsSet.add_agent(agent)
                            if not agent._everhighrisk_bool:
                                self.NewHRrolls.add_agent(agent)

                            agent._mean_num_partners = (
                                agent._mean_num_partners + params.HR_partnerScale
                            )
                            agent._highrisk_bool = True
                            agent._everhighrisk_bool = True
                            agent._highrisk_time = params.HR_M_dur

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

        elif (
            self.runRandom.random()
            < params.DemographicParams[race_type][sex_type]["INCAR"]
            * (1 + (hiv_bool * 4))
            * params.cal_IncarP
        ):
            if agent._SO == "HF":
                jailDuration = prob.HF_jail_duration
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
            self.totalIncarcerated += 1

            # PUT PARTNERS IN HIGH RISK
            for tmpA in agent._partners:
                if tmpA._highrisk_bool is True:
                    pass
                else:
                    if self.runRandom.random() < params.HR_proportion:
                        if not tmpA._highrisk_bool:
                            self.highrisk_agentsSet.add_agent(tmpA)
                            if not tmpA._everhighrisk_bool:
                                self.NewHRrolls.add_agent(tmpA)
                            tmpA._mean_num_partners += (
                                params.HR_partnerScale
                            )  # 32.5 #2 + 3.25 from incar HR
                            tmpA._highrisk_bool = True
                            tmpA._everhighrisk_bool = True
                            tmpA._highrisk_time = params.HR_F_dur
                if params.flag_PrEP and (
                    params.PrEP_target_model == "Incar"
                    or params.PrEP_target_model == "IncarHR"
                ):
                    # Atempt to put partner on prep if less than probability
                    if not tmpA._HIV_bool and agent.vaccine_time == 0:
                        self._initiate_PrEP(tmpA, time)

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

        if not tested:
            test_prob = params.DemographicParams[race_type][sex_type]["HIVTEST"]

            # Rescale based on calibration param
            test_prob = test_prob * params.cal_TestFreq

            # If roll less than test probablity
            if self.runRandom.random() < test_prob:
                # Become tested, add to tested agent set
                agent._tested = True
                self.NewDiagnosis.add_agent(agent)
                self.Trt_Tstd_agentSet.add_agent(agent)
                # If treatment co-enrollment enabled and coverage greater than 0
                if self.treatmentEnrolled and params.treatmentCov > 0: #TODO fix this logic
                    # For each partner, attempt to test for HIV
                    for ptnr in agent._partners:
                        if ptnr._HIV_bool and not ptnr._tested:
                            if self.runRandom.random() < 0.87:
                                ptnr._tested = True
                                self.NewDiagnosis.add_agent(ptnr)

    def _initiate_HAART(self, agent: Agent, time: int):
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
        assert agent._HIV_bool is True

        agent_haart = agent._HAART_bool
        agent_Test_bool = agent._tested
        agent_race = agent._race
        agent_so = agent._SO

        # Set HAART agents adherence at t=0 for all instanced HAART
        if time == 0 and agent_haart:
            agent_haart_adh = agent._HAART_adh
            if agent_haart_adh == 0:
                tmp_rnd = self.runRandom.random()
                HAART_ADH = params.DemographicParams[agent_race][agent_so]["HAARTadh"]
                if tmp_rnd < HAART_ADH:
                    adherence = 5
                else:
                    adherence = self.runRandom.randint(1, 4)

                # add to agent haart set
                agent._HAART_adh = adherence
                agent._HAART_time = time

        # Determine probability of HIV treatment
        if time >= 0 and agent_Test_bool:
            # Go on HAART
            if not agent_haart and agent._HAART_time == 0:
                if (
                    self.runRandom.random()
                    < params.DemographicParams[agent_race][agent_so]["HAARTprev"]
                    * params.cal_ART_cov
                ):

                    tmp_rnd = self.runRandom.random()
                    HAART_ADH = params.DemographicParams[agent_race][agent_so][
                        "HAARTadh"
                    ]
                    if tmp_rnd < HAART_ADH:
                        adherence = 5
                    else:
                        adherence = self.runRandom.randint(1, 4)

                    # Add agent to HAART class set, update agent params
                    agent._HAART_bool = True
                    agent._HAART_adh = adherence
                    agent._HAART_time = time
                    self.Trt_ART_agentSet.add_agent(agent)
            elif (
                agent_haart
                and self.runRandom.random()
                < params.DemographicParams[agent_race][agent_so]["HAARTdisc"]
            ):
                if agent._incar_treatment_time > 0 and params.inc_treat_RIC:
                    pass
                else:
                    agent._HAART_bool = False
                    agent._HAART_adh = 0
                    agent._HAART_time = 0
                    self.Trt_ART_agentSet.remove_agent(agent)

    def _PrEP_eligible(
        self, agent: Agent, time: int
    ) -> bool:  # REVIEWED should this be in agent? self not used - move to agent
        eligible = False
        if params.PrEP_target_model == "Allcomers":
            eligible = True
        elif params.PrEP_target_model == "CDCwomen":
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
        elif params.PrEP_target_model == "CDCmsm":
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
        elif params.PrEP_target_model == "HighPN10":
            if agent._mean_num_partners >= 10:
                eligible = True
        elif params.PrEP_target_model == "SRIns":
            if agent._sexualRole == "Insertive":
                eligible = True
        elif params.PrEP_target_model == "MSM":
            if agent._SO == ("MSM" or "MTF"):
                eligible = True
        elif params.PrEP_target_model == "RandomTrial":
            # If using random trial
            if time == 0:
                # if in init timestep 0, use agent set elligiblity
                eligible = (
                    agent._PrEP_eligible
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
        # If force flag set, auto kick off prep.
        if force is True:
            self.Trt_PrEP_agentSet.remove_agent(agent)
            agent._PrEP_bool = False
            agent._PrEP_reason = []
        # else if agent is no longer enrolled on PrEP, increase time since last dose
        elif agent._PrEP_time > 0:
            if agent._PrEP_time == 1:
                agent._PrEP_bool = False
                agent._PrEP_reason = []
                agent._PrEP_time -= 1
            else:
                agent._PrEP_time -= 1

        # else if agent is on PrEP, see if they should discontinue
        elif agent._PrEP_bool and agent._PrEP_time == 0:
            if (
                self.runRandom.random()
                < params.DemographicParams[agent._race][agent._SO]["PrEPdisc"]
            ):
                agent._PrEP_time = params.PrEP_falloutT
                self.Trt_PrEP_agentSet.remove_agent(agent)

                if "Oral" in params.PrEP_type:
                    agent._PrEP_bool = False
                    agent._PrEP_reason = []
            else:  # if not discontinue, see if its time for a new shot.
                if agent._PrEP_lastDose > 2:
                    agent._PrEP_lastDose = -1

        if params.PrEP_type == "Inj":
            self._calc_PrEP_load(agent)

    def initiate_Vaccine(self, agent: Agent, time: int, vaxType: str):
        """
        :Purpose:
            Progress vaccine. Agents may receive injection or progress in time since injection.

        :Input:
            agent: Agent
            time: int

        :Output:
            none
        """
        timeSinceVaccination = agent.vaccine_time - 1

        def vaccinate(ag: Agent, vax: str):
            ag.vaccine_time = 1
            ag.vaccine_type = vax


        if time == 1:
            if self.runRandom.random() < params.DemographicParams[agent._race][agent._SO]["vaccinePrev"]:
                vaccinate(agent, vaxType)
        else:
            if agent.vaccine_time >= 1:
                agent.vaccine_time += 1
            if (
                params.flag_booster
                and timeSinceVaccination
                % params.DemographicParams[agent._race][agent._SO]["boosterInterval"]
                == 1
                and  self.runRandom.random() < params.DemographicParams[agent._race][agent._SO]["boosterProb"]
            ):
                vaccinate(agent, vaxType)

        if time >= params.vaccine_start:
            if (
                random.random()
                < params.DemographicParams[agent._race][agent._SO]["prevVaccine"]
                * params.cal_Vaccine
            ):
                vaccinate(agent, vaxType)

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

        def _enrollPrEP(self, agent: Agent):
            agent._PrEP_bool = True
            agent._PrEP_time = 0
            self.Trt_PrEP_agentSet.add_agent(agent)
            self.newPrEPagents.add_agent(agent)
            self.newPrEPenrolls += 1
            if params.PrEP_target_model == "CDCwomen":
                if "IDU" in agent._PrEP_reason:
                    self.IDUprep += 1
                if "HIV test" in agent._PrEP_reason:
                    self.HIVprep += 1
                if "MSMW" in agent._PrEP_reason:
                    self.MSMWprep += 1
            tmp_rnd = self.runRandom.random()
            if params.PrEP_Adherence == "AtlantaMSM":
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
            if params.PrEP_type == "Inj":
                agent._PrEP_load = params.PrEP_peakLoad
                agent._PrEP_lastDose = 0

        # agent must exist
        assert agent is not None

        # Check valid input
        # Prep only valid for agents not on prep and are HIV negative
        if agent._PrEP_bool or agent._HIV_bool or agent.vaccine_time == 0:
            return None

        # Determine probability of HIV treatment
        agent_race = agent._race

        if force:
            _enrollPrEP(self, agent)
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
                numPrEP_agents < target_PrEP
                and time >= params.PrEP_startT
                and self._PrEP_eligible(agent, time)
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

        iterations = 1
        while iterations < 3:
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
                        % (minNum, maxNum, iterations)
                    )
                )
                iterations += 1

        print("No suitable PrEP agent")
        return None

    def _progress_to_AIDS(self, agent: Agent, agent_drug_type: str):
        """
        :Purpose:
            Model the progression of HIV agents to AIDS agents
        """
        # only valid for HIV agents
        if not agent._HIV_bool:
            raise ValueError(
                "HAART only valid for HIV agents!agent:%s" % str(agent._ID)
            )

        # if agent not in self.AIDS_agents:
        if not agent._HAART_bool:
            adherenceStat = agent._HAART_adh
            p = prob.adherence_prob(adherenceStat)

            if self.runRandom.random() < p * params.cal_ProgAIDS:
                agent._AIDS_bool = True
                self.HIV_AIDS_agentSet.add_agent(agent)

    def _reset_death_count(self):
        self.num_Deaths = {}
        self.deathSet = []
        for HIV_status in ["Total", "HIV-", "HIV+"]:
            self.num_Deaths.update({HIV_status: {}})
            for tmp_type in [HIV_status, "MSM", "HM", "HF", "WSW", "MTF"]:
                self.num_Deaths[HIV_status].update({tmp_type: 0})

    def _die_and_replace(self, time: int):
        """
        :Purpose:
            Let agents die and replace the dead agent with a new agent randomly.
        """
        totalDeaths = 0
        for agent in self.All_agentSet._members:

            if agent._incar_bool:
                pass
            else:
                # Probability for dying
                sex_type = agent._SO
                HIV_status = agent._HIV_bool
                AIDSStatus = agent._AIDS_bool
                agent_Race = agent._race

                if HIV_status:
                    if AIDSStatus:  # AIDS DEATH RATE
                        if agent_Race == "WHITE":
                            p = 34.4
                        elif agent_Race == "BLACK":
                            p = 41.6
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))

                    elif agent._HAART_adh > 1:  # HAART DEATH RATE
                        if agent_Race == "WHITE":
                            p = 8.6
                        elif agent_Race == "BLACK":
                            p = 10.4
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))

                    else:  # HIV+ DEATH RATE
                        if agent_Race == "WHITE":
                            p = 17.2
                        elif agent_Race == "BLACK":
                            p = 20.8
                        else:
                            raise ValueError("Invalid RACE type! %s" % str(agent_Race))
                    p = p * params.cal_Mortality

                elif not HIV_status:  # NON HIV DEATH RATE
                    if agent_Race == "WHITE":
                        p = 8.6
                    elif agent_Race == "BLACK":
                        p = 10.4
                    else:
                        raise ValueError("Invalid RACE type! %s" % str(agent_Race))

                else:
                    raise ValueError("Invalid HIV type! %s" % str(HIV_status))

                p = (
                    p / 12000.0
                )  # putting it into per 1 person-month from per 1000 person years

                if self.runRandom.random() < p:
                    totalDeaths += 1
                    if HIV_status:
                        ident = "HIV+"
                    else:
                        ident = "HIV-"

                    self.num_Deaths["Total"][sex_type] += 1
                    self.num_Deaths[ident][sex_type] += 1
                    self.deathSet.append(agent)
                    ID_number = agent.get_ID()
                    race = agent._race

                    # End all existing relationships
                    for rel in agent._relationships:
                        rel.progress(forceKill=True)

                        self.Relationships.remove_agent(rel)

                    # Remove agent node and edges from network graph
                    self.get_Graph().remove_node(agent)

                    # Remove agent from agent class and sub-sets
                    self.All_agentSet.remove_agent(agent)

                    # Delete agent object
                    del agent

                    # Create new agent
                    agent_cl = self._return_new_Agent_class(ID_number, race, sex_type)
                    self.create_agent(agent_cl, race)
                    self.get_Graph().add_node(agent_cl)

    def return_results(self):
        return self.ResultDict
