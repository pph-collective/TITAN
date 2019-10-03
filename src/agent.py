#!/usr/bin/env python
# encoding: utf-8

import random


class Agent:
    "Class for agent objects."

    def __init__(self, ID, SO, age, race, DU, initial_agent=False):
        """
        :Purpose:
            Constructor for an Agent

        :Input:
            :ID: Unique ID number to track each person agent
            :SO: Sexual Orientation of the agent
            :age: Age of the agent
            :race: Race of the agent (currently just "BLACK" or "WHITE")
            :DU: Drug use type ("IDU", "NIDU", "NDU")
            :initial_agent: whether this is the first agent created (default is False)
        """

        # self._ID is unique ID number used to track each person agent.
        self._ID = ID
        self._timeAlive = 0

        # self._initial_agent is set to "True" for agents that were used to
        # initialize the model.
        self._initial_agent = initial_agent

        # _parent_agent stores the parent agent if this agent is a member of an
        # Agent_set class instance.
        self._parent_agent = None

        # agent properties #REVIEW the `_` implies these are private, but they're used all over
        self._SO = SO
        self._age = age
        self._race = race
        self._DU = DU
        self._gender = SO[
            -1:
        ]  # Takes last letter of HM, HF, MSM, WSW, BiM, BiF to get agent gender. CANT USE MSMW! #REVIEW - inconsistent female

        self._ageBin = 0
        self._MSMW = False

        # agent-partner params
        self._relationships = []
        self._partners = []
        self._num_sex_partners = 0
        self._num_NE_partners = 0
        self._mean_num_partners = 0
        self._sexualRole = "Vers"

        # agent STI params
        self._STI_pool = []
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
        self._OAT_bool = False
        self._naltrex_bool = False
        self._DOC_OAT_bool = False
        self._DOC_NAL_bool = False
        self._MATprev = 0
        self._oatValue = 0
        self._PrEP_reason = []

        # PrEP pharmacokinetics
        self._PrEP_load = 0.0
        self._PrEP_lastDose = 0

        # agent high risk params
        self._highrisk_type = None
        self._highrisk_bool = False
        self._highrisk_time = 0
        self._everhighrisk_bool = False

        # agent incarcartion params
        self._incar_bool = False
        self._ever_incar_bool = False
        self._incar_time = 0
        self._incar_treatment_time = 0

    def __str__(self):
        return "\t%.6d\t%d\t%s\t%s\t%s\t%s\t%s" % (
            self._ID,
            self._age,
            self._gender,
            self._SO,
            self._DU,
            self._race,
            self._HIV_bool,
        )

    def __repr__(self):
        return str(self._ID)

    def get_ID(self):
        return self._ID

    def set_parent_agent(self, agent):
        self._parent_agent = agent

    def get_parent_agent(self):
        return self._parent_agent

    def bond(self, partner, relationship=None):
        if relationship == None:
            exit(9)
        self._relationships.append(relationship)
        partner._relationships.append(relationship)

        # Test this new code!!!
        self.pair(partner)
        partner.pair(self)

        self._num_sex_partners += 1
        partner._num_sex_partners += 1

    def unbond(self, partner, relationship=None):
        if relationship == None:
            exit(9)
        self._relationships.remove(relationship)
        partner._relationships.remove(relationship)

        self.unpair(relationship._ID1)
        self.unpair(relationship._ID2)
        partner.unpair(relationship._ID1)
        partner.unpair(relationship._ID2)

    def pair(self, partner):
        "Pairs two agents to each other."
        if partner.get_ID() != self.get_ID():
            if partner in self._partners:
                print("ASDF")
            else:
                self._partners.append(partner)

    def unpair(self, agent):
        "Removes agent from agent set."
        try:
            if self != agent:
                self._partners.remove(agent)
        except KeyError:
            raise KeyError(
                "agent %s is not a member of agent set %s"
                % (agent.get_ID(), self.get_ID())
            )

    def partner_list(self):
        ptnrs = list()
        if self._partners != None:
            for temp in self._partners:
                ptnrs.append(temp._ID)

        return ptnrs

    # REVIEW seems like this should reference __str__ logic instead of rewriting it
    def print_agent(self):
        print(
            "\t%.6d\t%d\t%s\t%s\t%s\t%s\t%s\t%s"
            % (
                self._ID,
                self._age,
                self._gender,
                self._SO,
                self._DU,
                self._race,
                self._HIV_bool,
                self._incar_bool,
            )
        )

    def print_agent_abridge(self):
        print("\t%.6d\t%s\t%s\t%s" % (self._ID, self._gender, self._SO, self._DU))

    def vars(self):
        return "%.6d,%d,%s,%s,%s,%s,%s,%d,%d,%d,%s,%s,%s,%s\n" % (
            self._ID,
            self._age,
            self._gender,
            self._SO,
            self._DU,
            self._race,
            self._HIV_bool,
            len(self._partners),
            self._num_sex_partners,
            self._timeAlive,
            self._AIDS_bool,
            self._tested,
            self._PrEP_bool,
            self._incar_bool,
        )

    def print_agent_to_file(self, filename, time=None, overWrite="a"):
        if overWrite == "a":
            agentList = ""
        else:
            agentList = "Time,ID,Age,Gdr,SO,DU,Race,HIV+,Ptnrs,AIDS,Tested,PrEP,Incar\n"

        agentList += str(time) + "," + self.vars()

        open(str(filename), overWrite).write(agentList)

    def print_relationships(self):
        for tmpR in self._relationships:
            tmpR.print_rel()


class Relationship:
    "Class for agent relationships."

    def __init__(self, ID1, ID2, SO, rel_type, duration, initial_agent=False):
        """
        :Purpose:
            Constructor for a Relationship

        :Input:
            :ID1: ID of first agent
            :ID2: ID of second agent
            :SO: Orientation of relationship
            :rel_type: #REVIEW
            :duration: length of relationship
            :initial_agent: #REVIEW
        """

        # self._ID is unique ID number used to track each person agent.
        self._ID1 = ID1
        self._ID2 = ID2
        self._ID = self._ID1.get_ID() * 100000 + self._ID2.get_ID()
        # self._initial_agent is set to "True" for agents that were used to
        # initialize the model.
        self._initial_agent = initial_agent

        # Relationship properties
        self._SO = SO
        self._rel_type = rel_type
        self._duration = duration
        self._total_sex_acts = 0

        # Relationship STI params
        self._STI_pool = []
        self._HIV_bool = False
        self._HIV_time = 0
        self._AIDS_bool = False
        self._AIDS_time = 0

        # Relationship treatment params
        self._tested = False
        self._HAART_bool = False
        self._HAART_time = 0
        self._HAART_adh = 0

        # Relationship incarcartion params
        self._incar_bool = False
        self._incar_time = 0

    def progress(self, forceKill=False):
        if self._duration <= 0 or forceKill:
            agent = self._ID1
            partner = self._ID2

            agent.unbond(partner, self)

            return True
        else:
            self._duration = self._duration - 1
            return False

    def get_ID(self):
        return self._IDq

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


class Agent_set(Agent):
    """
    Class for agents that contain a "set" of agents from a lower
    hierarchical  level.
    """

    def __init__(self, world, ID, parent=None, numerator=None):
        # _members stores agent set members in a dictionary keyed by ID
        self._ID = ID
        self._members = []
        self._subset = {}

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

    def clear_set(self):
        self._members = []
        self._subset = {}

    def get_agents(self):
        return self._members.__iter__()

    def get_agent(self, ID):
        "Returns an agent given the agent's ID"
        return self._members[ID]

    def is_member(self, ID):
        "Returns true if agent is a member of this set"
        return ID in self._members

    def add_agent(self, agent):
        "Adds a new agent to the set."
        self._members.append(agent)

    def remove_agent(self, agent):
        "Removes agent from agent set."
        try:
            self._members.remove(agent)
        except:
            pass

        for tmpS in self.iter_subset():
            tmpS.remove_agent(agent)

    def iter_agents(self):
        for agent in self.get_agents():
            yield agent

    def num_members(self):
        return len(self._members)

    def random_agent(self):
        try:
            return random.choice(self._members)
        except:
            pass

    def add_subset(self, subset):
        "Adds a new Agent_set to the current sets subset."
        if subset._ID in self._subset:
            raise KeyError(
                "Subset %s is already a member of Agent_set set %s"
                % (subset.get_ID(), self._ID)
            )
        self._subset[subset._ID] = subset

        # Set the agent's _parent_agent to reflect the parent of this Agent_set
        # instance (self)
        subset.set_parent_agent(self)

    def remove_subset(self, subset):
        "Removes Agent_set to the current sets subset."
        try:
            self._subset.pop(subset.ID)
        except KeyError:
            raise KeyError(
                "subset %s is not a member of set %s" % (subset.get_ID(), self.get_ID())
            )

    def iter_subset(self):
        for subset in list(self._subset.values()):
            yield subset

    def set_parent_set(self, master_set):
        self._parent_set = master_set

    def get_parent_set(self):
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
