#!/usr/bin/env python
# encoding: utf-8

from typing import Sequence, List, Dict, Optional
from . import params


class Agent:
    "Class for agent objects."

    # class variable for agent creation
    next_agent_id = 0

    @classmethod
    def update_id_counter(cls):
        cls.next_agent_id += 1

    def __init__(
        self, SO: str, age: int, race: str, DU: str, initial_agent: bool = False,
    ):
        """
        Initialize an agent based on given properties

        args:
            ID (int) - Unique agent ID
            SO (str) - Sexual orientation flag (HM, HF, MSM)
            age (int) - Agents initialization age
            race (str) - Race of agent
            DU (str) - Drug use flag (IDU, NIDU, NDU)
            initial_agent (bool) - If the agent was created during model init

        returns:
            None
        """
        # self._ID is unique ID number used to track each person agent.
        self._ID = self.next_agent_id
        self.update_id_counter()
        self._timeAlive = 0

        # agent properties
        self._SO = SO
        self._age = age
        self._race = race
        self._DU = DU

        self._ageBin = 0
        self._MSMW = False

        # agent-partner params
        self._relationships: List[Relationship] = []
        self._partners: List[Agent] = []
        self._num_sex_partners = 0  # TO_REVIEW - not actually used anywhere
        self._num_NE_partners = 0  # TO_REVIEW - not actually used anywhere
        self._mean_num_partners = 0
        self._sexualRole = "Vers"

        # agent STI params
        self._HIV_bool = False
        self._HIV_time = 0
        self._AIDS_bool = False
        self._AIDS_time = 0
        self._PrEPresistance = 0

        # agent treatment params
        self._tested = False
        self._HAART_bool = False
        self._HAART_time = 0
        self._HAART_adh = 0
        self._SNE_bool = False
        self._PrEP_bool = False
        self._PrEP_time = 0
        self._PrEP_adh = 0
        self._treatment_bool = False
        self._treatment_time = 0
        self._PrEP_reason: List[str] = []

        # PrEP pharmacokinetics
        self._PrEP_load = 0.0
        self._PrEP_lastDose = 0

        # agent high risk params
        self._highrisk_bool = False
        self._highrisk_time = 0
        self._everhighrisk_bool = False

        # agent incarcartion params
        self._incar_bool = False
        self._ever_incar_bool = False
        self._incar_time = 0
        self._incar_treatment_time = 0

    def __str__(self):
        """
        String formatting of agent object

        returns:
            String formatted tab-deliminated agent properties
        """
        return "\t%.6d\t%d\t%s\t%s\t%s\t%s" % (
            self._ID,
            self._age,
            self._SO,
            self._DU,
            self._race,
            self._HIV_bool,
        )

    def __repr__(self):
        """
        Repr formatting of agent object

        returns:
            ID (str) - agent ID as str
        """
        return str(self._ID)

    def get_ID(self):
        """
        Get the agent ID

        returns:
            ID (int) - agent ID
        """
        return self._ID

    def partner_list(self):
        """
        Return the list of partners for this agent

        returns:
            _partners (list) - list of partners
        """
        ptnrs = list()
        if self._partners is not None:
            for temp in self._partners:
                ptnrs.append(temp._ID)

        return ptnrs

    def get_acute_status(self, time: int) -> bool:
        """
        :Purpose:
            Get acute status of agent at time period
        :Input:
            time : int
        :Output:
            acute_status : bool
        """
        acuteTimePeriod = 2
        hiv_t = self._HIV_time

        if hiv_t <= acuteTimePeriod and hiv_t > 0:
            return True
        else:
            return False

    def PrEP_eligible(self, time: int) -> bool:
        """
        :Purpose:
            Determine if an agent is eligible for PrEP

        :Input:
            time : int

        :Output:
            eligibility : bool
        """
        eligible = False
        if params.PrEP_target_model == "Allcomers":
            eligible = True
        elif params.PrEP_target_model == "CDCwomen":
            if self._SO == "HF":
                for rel in set(self._relationships):
                    if rel._ID1 == self:
                        partner = rel._ID2
                    else:
                        partner = rel._ID1

                    if rel._duration > 1:
                        if partner._DU == "IDU":
                            eligible = True
                            self._PrEP_reason.append("IDU")
                        if partner._tested:
                            eligible = True
                            self._PrEP_reason.append("HIV test")
                        if partner._MSMW:
                            eligible = True
                            self._PrEP_reason.append("MSMW")
        elif params.PrEP_target_model == "CDCmsm":
            if self._SO == "MSM":
                for rel in self._relationships:
                    if rel._ID1 == self:
                        partner = rel._ID2
                    else:
                        partner = rel._ID1

                    if rel._duration > 1:
                        if partner._tested or self._mean_num_partners > 1:
                            eligible = True
        elif params.PrEP_target_model == "HighPN5":
            if self._mean_num_partners >= 5:
                eligible = True
        elif params.PrEP_target_model == "HighPN10":
            if self._mean_num_partners >= 10:
                eligible = True
        elif params.PrEP_target_model == "SRIns":
            if self._sexualRole == "Insertive":
                eligible = True
        elif params.PrEP_target_model == "MSM":
            if self._SO in ("MSM", "MTF"):
                eligible = True

        return eligible

    def update_PrEP_load(self):
        """
        :Purpose:
            Determine and update load of PrEP concentration in agent.

        :Input:
            none

        :Output:
            none
        """
        # N(t) = N0 (0.5)^(t/t_half)
        self._PrEP_lastDose += 1
        if self._PrEP_lastDose > 12:
            self._PrEP_load = 0.0
        else:
            self._PrEP_load = params.PrEP_peakLoad * (
                (0.5) ** (self._PrEP_lastDose / (params.PrEP_halflife / 30))
            )

    def get_transmission_probability(self, interaction: str) -> float:
        """ Decriptor
        :Purpose:
            Determines the probability of a transmission event based on interaction type.

        :Input:
            interaction : str - "NEEDLE" or "SEX"

        :Output:
            probability : float
        """

        sex_type = self._SO
        agentAdherence = self._HAART_adh

        # Logic for if needle or sex type interaction
        p: float
        if interaction == "NEEDLE":
            p = params.TransmissionProbabilities["NEEDLE"][str(agentAdherence)]
        elif interaction == "SEX":
            p = params.TransmissionProbabilities["SEX"][sex_type][str(agentAdherence)]

        # Scaling parameter for acute HIV infections
        if self.get_acute_status(0):
            p = p * params.cal_AcuteScaling

        # Scaling parameter for positively identified HIV agents
        if self._tested:
            p = p * (1 - params.cal_RR_Dx)

        # Tuning parameter for ART efficiency
        if self._HAART_bool:  # self.AdherenceAgents[agent] > 0:
            p = p * params.cal_RR_HAART

        # Racial calibration parameter to attain proper race incidence disparity
        if self._race == "BLACK":
            p = p * params.cal_raceXmission

        # Scaling parameter for per act transmission.
        p = p * params.cal_pXmissionScaling

        return p

    # REVIEWED as written, it's not guaranteed to return an int, but seems likely given params, can't this be more clearly written as a for loop? it also looks like the params adds the p values, while this loop also adds p-values, is that correct? - Sarah to review with Max
    def get_number_of_sexActs(self, rand_gen) -> int:
        """
        :Purpose:
            Number of sexActs an agent has done.

        :Input:
            rand_gen : random number generator (e.g. self.runRandom in ABM_core)

        :Output:
            number_sex_act : int
        """
        # 1 time per year 96 1.9 29 0.9 67 3.4
        # 2–5 times per year 428 8.2 184 5.8 244 12.2
        # 6–11 times per year 328 6.3 183 5.7 145 7.3
        # 12–23 times per year 376 7.2 251 7.9 125 6.3
        # 24–35 times per year 1,551 29.9 648 20.3 903 45.3
        # 36–51 times per year 1,037 20.0 668 20.9 369 18.5
        # 52–155 times per year 644 12.4 605 18.9 39 2.0
        # >156 times per year 733 14.1 631 19.7 102 5.1
        rv = rand_gen.random()
        pMatch = 0.0

        for i in range(1, 6):
            pMatch = params.sexualFrequency[i]["p_value"]
            if rv <= pMatch:
                minSA = params.sexualFrequency[i]["min"]
                maxSA = params.sexualFrequency[i]["max"]
                return rand_gen.randrange(minSA, maxSA, 1)

        # REVIEWED - fallthrough return? 0? - Sarah to review with Max
        return 0


class Relationship:
    "Class for agent relationships."

    # class variable for relationship creation
    next_rel_id = 0

    @classmethod
    def update_id_counter(cls):
        cls.next_rel_id += 1

    def __init__(self, ID1: Agent, ID2: Agent, SO: str, rel_type: str, duration: int):
        """
        :Purpose:
            Constructor for a Relationship

        :Input:
            :ID1: first agent
            :ID2: second agent
            :SO: Orientation of relationship
            :rel_type: (future feature) - sex bond or idu bond or both
            :duration: length of relationship
        """
        # make sure these agents can be in a relationship
        assert ID1 != ID2, "Cannot create relationship with same agent"
        assert (
            ID1 not in ID2._partners and ID2 not in ID1._partners
        ), "Agent's already partnered!"

        # self._ID is unique ID number used to track each person agent.
        self._ID1 = ID1
        self._ID2 = ID2
        self._ID = self.next_rel_id
        self.update_id_counter()

        # Relationship properties
        self._duration = duration
        self._total_sex_acts = 0

        self.bond(ID1, ID2)

    def progress(self, forceKill: bool = False):
        if self._duration <= 0 or forceKill:
            self.unbond(self._ID1, self._ID2)
            return True
        else:
            self._duration = self._duration - 1
            return False

    def bond(self, agent: "Agent", partner: "Agent"):
        """
        Bond two agents. Adds each to a relationship object, then partners in each
        others' partner list.

        args:
            agent: Agent - new partner of partner
            partner: Agent - new partner of agent
        returns:
            None
        """

        # Append relationship to relationships list for each agent
        agent._relationships.append(self)
        partner._relationships.append(self)

        # Pair agent with partner and partner with agent
        agent._partners.append(partner)
        partner._partners.append(agent)

        # Increment number of sex partners for both
        agent._num_sex_partners += 1
        partner._num_sex_partners += 1

    def unbond(self, agent: "Agent", partner: "Agent"):
        """
        Unbond two agents. Removes relationship from relationship lists. Removes partners in each others' partner list.

        args:
            agent: Agent - former partner of partner
            partner: Agent - former partner of agent
        returns:
            None
        """

        # Remove relationship to relationships list for each agent
        agent._relationships.remove(self)
        partner._relationships.remove(self)

        # Unpair agent with partner and partner with agent
        agent._partners.remove(partner)
        partner._partners.remove(agent)

    def get_ID(self):
        return self._ID

    def __repr__(self):
        return "\t%.6d\t%.6d\t%s\t%s\t%d\t%d" % (
            self._ID1.get_ID(),
            self._ID2.get_ID(),
            self._SO,
            self._rel_type,
            self._duration,
            self._total_sex_acts,
        )

    def print_rel(self):
        return "\t%.6d\t%.6d\t%s\t%s\t%d\t%d" % (
            self._ID1.get_ID(),
            self._ID2.get_ID(),
            self._SO,
            self._rel_type,
            self._duration,
            self._total_sex_acts,
        )

    def print_rels(self):
        print("\t_____________ %s _____________" % self.get_ID())
        print("\tID1\tID2\tSO\tRel\tDur\tSexA")

        for tmpA in self.iter_agents():
            print(tmpA)
        print("\t______________ END ______________")

    def vars(self):
        return "%.6d,%.6d,%s,%s,%d,%d\n" % (
            self._ID1.get_ID(),
            self._ID2.get_ID(),
            self._SO,
            self._rel_type,
            self._duration,
            self._total_sex_acts,
        )


class Agent_set:
    """
    Class for agents that contain a "set" of agents from a lower
    hierarchical  level.
    """

    def __init__(
        self, ID: str, parent: "Agent_set" = None, numerator: "Agent_set" = None
    ):
        # _members stores agent set members in a dictionary keyed by ID
        self._ID = ID
        self._members: List[Agent] = []
        self._subset: Dict[str, Agent_set] = {}

        # _parent_set stores the parent set if this set is a member of an
        # Agent_set class instance. For example, for a set that is a
        # member of a larger set, the _parent_set for that set  would
        # be that larger set.
        self._parent_set = parent
        if parent:
            parent.add_subset(self)
        if numerator:
            self._numerator = numerator
        else:
            self._numerator = self

    def __repr__(self):
        return self._ID

    def __str__(self):
        return self._ID

    def clear_set(self):
        self._members: List[Agent] = []
        self._subset: Dict[str, str] = {}

    def get_ID(self):
        return self._ID

    def get_agents(self):
        return self._members.__iter__()

    def is_member(self, agent: Agent):
        "Returns true if agent is a member of this set"
        return agent in self._members

    def add_agent(self, agent: Agent):
        "Adds a new agent to the set."
        self._members.append(agent)

    def remove_agent(self, agent: Agent):
        "Removes agent from agent set."
        if self.is_member(agent):
            self._members.remove(agent)

        for tmpS in self.iter_subset():
            tmpS.remove_agent(agent)

    def iter_agents(
        self,
    ):  # REVIEWED isn't this redundant with get_agents? why not have __iter__ return the agents? then we could use the syntax agent in agent_set - maybe consolidate later
        for agent in self.get_agents():
            yield agent

    def num_members(self) -> int:
        return len(self._members)

    def add_subset(self, subset: "Agent_set"):
        "Adds a new Agent_set to the current sets subset."
        if subset._ID in self._subset:
            raise KeyError(
                "Subset %s is already a member of Agent_set set %s"
                % (subset.get_ID(), self._ID)
            )
        self._subset[subset._ID] = subset

    def remove_subset(self, subset: "Agent_set"):
        "Removes Agent_set to the current sets subset."
        try:
            self._subset.pop(subset._ID)
        except KeyError:
            raise KeyError(
                "subset %s is not a member of set %s" % (subset.get_ID(), self.get_ID())
            )

    def iter_subset(self):
        for subset in list(self._subset.values()):
            yield subset

    def set_parent_set(self, master_set: "Agent_set"):
        self._parent_set = master_set

    def get_parent_set(self) -> Optional["Agent_set"]:
        return self._parent_set

    def print_subsets(self):
        print("\t__________ %s __________" % self.get_ID())
        print("\tID\t\tN\t\t%")
        for tmpS in self.iter_subset():
            if tmpS.num_members() > 0:
                print(
                    "\t%-6s\t%-5d\t%.2f"
                    % (
                        tmpS._ID,
                        tmpS.num_members(),
                        (1.0 * tmpS.num_members() / tmpS._numerator.num_members()),
                    )
                )
            for tmpSS in tmpS.iter_subset():
                if tmpSS.num_members() > 0:
                    print(
                        "\t-%-4s\t%-5d\t%.2f"
                        % (
                            tmpSS._ID,
                            tmpSS.num_members(),
                            (
                                1.0
                                * tmpSS.num_members()
                                / tmpSS._numerator.num_members()
                            ),
                        )
                    )
        print("\t______________ END ______________")
