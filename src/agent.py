#!/usr/bin/env python
# encoding: utf-8
"""
*****************************************************************************
Author(s):	Maximilian King  (previous authors: Lars Seemann - lseemann@uh.edu)
Email: Maximilian_King@brown.edu
Organization: Marshall Lab, Department of Epidemiology - Brown University

Description:
    Contains classes to assist in making agents. 'Person' agents, for example,
would be subclasses of the Agent class, while Household, Neighborhood, and
Region agents would be represented as subclasses of the Agent_set object (as
households, neighborhoods, and regions all contain lower-level agents).


Copyright (c) 2016, Maximilian King
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************
"""


class Agent(object):
    "Class for agent objects."

    def __init__(self, ID, SO, age, race, DU, initial_agent=False):
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
        self._ID = ID
        self._timeAlive = 0

        # self._initial_agent is set to "True" for agents that were used to
        # initialize the model.
        self._initial_agent = initial_agent

        # _parent_agent stores the parent agent if this agent is a member of an
        # Agent_set class instance.
        self._parent_agent = None

        # agent properties
        self._SO = SO
        self._age = age
        self._race = race
        self._DU = DU
        self._gender = SO[
            -1:
        ]  # Takes last letter of HM, HF, MSM, WSW, BiM, BiF to get agent gender. CANT USE MSMW!

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

        # Set gender on small switch statement.

    def __str__(self):
        """
        String formatting of agent object

        returns:
            String formatted tab-deliminated agent properties
        """
        return "\t%.6d\t%d\t%s\t%s\t%s\t%s\t%s" % (
            self._ID,
            self._age,
            self._gender,
            self._SO,
            self._DU,
            self._race,
            self._HIV_bool,
        )
        # return str(self._ID)+ '\t'+str(self._SO)+'\t'+str(self._DU)+'\t'+ str(self._HIV_bool)

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

    def set_parent_agent(self, parent):
        """
        Set the parent of the agent.

        args:
            parent (Agent(object)) - Parent agent object
        returns:
            None
        """
        self._parent_agent = parent

    def get_parent_agent(self):
        """
        Get the parent of the agent.

        returns:
            parent (Agent(object)) - Parent agent object
        """
        return self._parent_agent

    def bond(self, partner, relationship=None):
        """
        Bond two agents. Adds each to a relationship object, then partners in each
        others' partner list.

        todo: Disentangle this from partner. These can be condensed.

        args:
            partner (Agent(object)) - new partner of agent
            relationship (Relationship(object)) - relationship for partnership
        returns:
            None
        """

        if relationship is None:
            exit(9)

        # Append relationship to relationships list for each agent
        self._relationships.append(relationship)
        partner._relationships.append(relationship)

        # Pair agent with partner and partner with agent
        self.pair(partner)
        partner.pair(self)

        # Increment number of sex partners for both
        self._num_sex_partners += 1
        partner._num_sex_partners += 1

    def unbond(self, partner, relationship=None):
        """
        Unbond two agents. Adds each to a relationship object, then partners in each
        others' partner list.

        todo: Disentangle this from partner. These can be condensed.

        args:
            partner (Agent(object)) - new partner of agent
            relationship (Relationship(object)) - relationship for partnership
        returns:
            None
        """
        if relationship is None:
            exit(9)

        # Remove relationship to relationships list for each agent
        self._relationships.remove(relationship)
        partner._relationships.remove(relationship)

        # Unpair agent with partner and partner with agent
        self.unpair(relationship._ID1)
        self.unpair(relationship._ID2)
        partner.unpair(relationship._ID1)
        partner.unpair(relationship._ID2)

    def pair(self, partner):
        """
        Pair two agents by adding each to the respective _partners list

        args:
            partner (Agent(object)) - new partner of agent
        returns:
            None
        """

        if partner.get_ID() != self.get_ID():
            if partner in self._partners:
                print(
                    "ASDF"
                )  # raise KeyError("Partner %s is already bonded with agent %s"%(partner.get_ID(), self._ID))
            # assert partner not in self._partners
            # todo: why was this (^) assertion and KeyError avoided?!
            else:
                self._partners.append(partner)

    def unpair(self, partner):
        """
        Unpair two agents by removing each to the respective _partners list

        args:
            agent (Agent(object)) - new partner of agent
        returns:
            None
        """
        try:
            if self != partner:
                self._partners.remove(partner)
        except KeyError:
            raise KeyError(
                "agent %s is not a member of agent set %s"
                % (partner.get_ID(), self.get_ID())
            )

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

    def print_agent(self):
        """
        Print the agent properties to stdout

        returns:
            None
        """
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
        """
        Print the abridged agent properties to stdout

        returns:
            None
        """
        print("\t%.6d\t%s\t%s\t%s" % (self._ID, self._gender, self._SO, self._DU))

    def vars(self):
        """
        Get agent specific variables (used for printing stats)

        returns:
            vars (str) - string of variables for agent
        """
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
        """
        Print the agent variables to a file

        args:
            filename (str) - name of file to write to
            time (int) - timestep of print
            overWrite (str) - write flag for f open
        returns:
            None
        """
        if overWrite == "a":
            agentList = ""
        else:
            agentList = "Time,ID,Age,Gdr,SO,DU,Race,HIV+,Ptnrs,AIDS,Tested,PrEP,Incar\n"

        agentList += str(time) + "," + self.vars()

        open(str(filename), overWrite).write(agentList)

    def print_relationships(self):
        """
        Print the agents relationships to stdout

        returns:
            None
        """
        for tmpR in self._relationships:
            tmpR.print_rel()


class Relationship(object):
    "Class for agent relationships."

    def __init__(self, ID1, ID2, SO, rel_type, duration, initial_agent=False):

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
            # print self._duration
            agent = self._ID1
            partner = self._ID2

            # print "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDeleting relationship between %d and %d" % (agent.get_ID(), partner.get_ID())
            # print "\tAgt:",agent
            # print "\tPtn:", partner
            # self.print_rel()
            agent.unbond(partner, self)
            # agent.unpair(partner)
            # partner.unpair(agent)
            # del self
            return True
            # self._ID1._relationships.remove(self)
            # self._ID2._relationships.remove(self)
        else:
            self._duration = self._duration - 1
            return False

    def delete(self):
        del self

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
        # print self.partner_list()
        print("\t_____________ %s _____________" % self.get_ID())
        print("\tID1\tID2\tSO\tRel\tDur\tSexA")

        for tmpA in self.iter_agents():
            print(tmpA)
        print("\t______________ END ______________")
        # print "\t%.6d\t%.6d\t%s\t%s\t%d\t%d" % (self._ID1.get_ID(), self._ID2.get_ID(), self._SO, self._rel_type, self._duration, self._total_sex_acts)

    def vars(self):
        return "%.6d,%.6d,%s,%s,%d,%d\n" % (
            self._ID1.get_ID(),
            self._ID2.get_ID(),
            self._SO,
            self._rel_type,
            self._duration,
            self._total_sex_acts,
        )


class Relationship_set(Relationship):
    """
    Class for agents that contain a "set" of relationsips from a lower
    hierarchical  level.
    """

    def __init__(self, world, ID, parent=None, numerator=None):
        # Agent.__init__(self, world, ID, initial_agent)

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


class Agent_set(Agent):
    """
    Class for agents that contain a "set" of agents from a lower
    hierarchical  level.
    """

    def __init__(self, world, ID, parent=None, numerator=None):
        # Agent.__init__(self, world, ID, initial_agent)

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
        # return self._members.values()
        return self._members.__iter__()

    def get_agent(self, ID):
        "Returns an agent given the agent's ID"
        return self._members[ID]

    def is_member(self, ID):
        "Returns true if agent is a member of this set"
        return ID in self._members

    def add_agent(self, agent):
        "Adds a new agent to the set."
        # agent.print_agent()
        # if agent in self._members:
        #    raise KeyError("agent %s is already a member of agent set %s"%(agent.get_ID(), self._ID))
        self._members.append(agent)
        # self._members[agent.get_ID()] = agent

        # if self._subset: #if subsets exist, try to add the agent from those sets
        #     try:
        #         self._subset[agent._SO].add_agent(agent)
        #     except:
        #         print "agent %s is already a member of agent set %s"%(agent.get_ID(), self.get_ID())
        # Set the agent's _parent_agent to reflect the parent of this Agent_set
        # instance (self)
        # agent.set_parent_agent(self)

    def remove_agent(self, agent):
        "Removes agent from agent set."

        # ID = agent.get_ID()
        # print "Removing agent %d"%ID

        try:
            self._members.remove(agent)
            # print "agent %s has been removed from agent set %s"%(agent._ID, self.get_ID())
            # self._members.
        except KeyError:
            # print "agent %s is not a member of agent set %s"%(agent.get_ID(), self.get_ID())
            # agent.print_agent()
            pass

        for tmpS in self.iter_subset():
            # tmpS.print_agents()

            tmpS.remove_agent(agent)
        # Reset the agent's _parent_agent
        # assert agent.get_parent_agent().get_ID() == self.get_ID(), "Removing agent from an Agent_set it does not appear to be assigned to."
        # agent.set_parent_agent(None)
        # print "Removed agent", agent._ID

    def iter_agents(self):
        for agent in self.get_agents():
            yield agent

    def num_members(self):
        return len(self._members)

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

    def print_agents(self):
        print("\t_____________ %s _____________" % self.get_ID())
        print("\tID\tAge\tGdr\tSO\tDU\tRace\tHIV+\tPtnrs")

        for tmpA in self.iter_agents():
            tmpA.print_agent()
        print("\t______________ END ______________")

    def print_agents_abridged(self):
        print("\t_________ %s _________" % self.get_ID())
        print("\tID\tGdr\tSO\tDU")

        for tmpA in self.iter_agents():
            tmpA.print_agent_abridge()
        print("\t__________ END __________")

    def print_agents_to_file(
        self, time=None, overWrite="a", filename="Results/Tableau_Agent_Output_File.txt"
    ):
        if overWrite == "a":
            agentList = ""
        else:
            agentList = "Time,ID,Age,Gdr,SO,DU,Race,HIV+,Ptnrs,LTPtnrs,AliveTime,AIDS,Tested,PrEP,Incar\n"
        for tmpA in self.iter_agents():
            agentList += str(time) + "," + tmpA.vars()

        open(filename, overWrite).write(agentList)

    def print_agent_relationshps(self):
        print("\t_____________ %s _____________" % self.get_ID())
        print("\tID1\t\tID2\t\tSO\tTy\tDur\tActs")

        for tmpR in self._members:
            tmpR.print_rel()
        print("\t______________ END ______________")

    def print_agent_relationships_to_file(self, time=None, overWrite="a"):
        if overWrite == "a":
            agentList = ""
        else:
            agentList = "Time,ID1,ID2,SO,Ty,Dur,Acts\n"
        for tmpR in self.iter_agents():
            agentList += str(time) + "," + tmpR.vars()

        open("Results/Tableau_Rel_Output_File.txt", overWrite).write(agentList)

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
