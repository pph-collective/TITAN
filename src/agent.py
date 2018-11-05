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

#import random

class Agent(object):
    "Class for agent objects."
    def __init__(self, ID, SO, age, race, DU, initial_agent=False):

        # self._ID is unique ID number used to track each person agent.
        self._ID = ID
        self._timeAlive = 0

        # self._initial_agent is set to "True" for agents that were used to
        # initialize the model.
        self._initial_agent = initial_agent

        # _parent_agent stores the parent agent if this agent is a member of an
        # Agent_set class instance. For example, for a person agent that is a
        # member of a household, the _parent_agent for that person agent would
        # be that household.
        self._parent_agent = None

        # agent properties
        self._SO = SO
        self._age = age
        self._race = race
        self._DU = DU
        self._gender = SO[-1:] #Takes last letter of HM, HF, MSM, WSW, BiM, BiF to get agent gender

        self._ageBin = 0

        # agent-partner params
        self._relationships = []
        self._partners = []
        self._num_sex_partners = 0
        self._num_NE_partners = 0
        self._mean_num_partners = 0
        self._sexualRole = 'Vers'

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

        #PrEP pharmacokinetics
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

        #Set gender on small switch statement.

    def __str__(self):
        return str(self._ID)+ '\t'+str(self._SO)+'\t'+str(self._DU)+'\t'+ str(self._HIV_bool)
    def __repr__(self):
        return str(self._ID)

    def get_ID(self):
        return self._ID

    def set_parent_agent(self, agent):
        self._parent_agent = agent

    def get_parent_agent(self):
        return self._parent_agent

    def bond(self, partner, relationship = None):
        #self._partners[0] = partner
        # self._partners.append(partner)
        # partner._partners.append(self)
        if relationship == None:
            exit(9)
        #tmp_relationship = Relationship(self, partner, "MSM", "SE", duration)
        self._relationships.append(relationship)
        partner._relationships.append(relationship)

        #Test this new code!!!
        self.pair(partner)
        partner.pair(self)

        self._num_sex_partners +=1
        partner._num_sex_partners +=1
        # self.pair(relationship._ID1)
        # self.pair(relationship._ID2)
        # partner.pair(relationship._ID1)
        # partner.pair(relationship._ID2)
        #partner.pair(self)
        #print "Agent %d bound to partner %d"%(self._ID, partner._ID)

    def unbond(self, partner, relationship = None):
        if relationship == None:
            exit(9)
        #tmp_relationship = Relationship(self, partner, "MSM", "SE", duration)
        #self.print_agent()
        #partner.print_agent()
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
                print "ASDF"#raise KeyError("Partner %s is already bonded with agent %s"%(partner.get_ID(), self._ID))
            #assert partner not in self._partners

            else:
                self._partners.append(partner)

        #print "Agent %d bound to partner %d"%(self._ID, partner._ID)
        # Set the agent's _parent_agent to reflect the parent of this Agent_set
        # instance (self)
        # agent.set_parent_agent(self)

    def unpair(self, agent):
        "Removes agent from agent set."
        try:
            if self != agent:self._partners.remove(agent)
            #self._members.
        except KeyError:
            raise KeyError("agent %s is not a member of agent set %s"%(agent.get_ID(), self.get_ID()))

    def partner_list(self):
        ptnrs = list()
        if self._partners != None:
            for temp in self._partners:
                ptnrs.append(temp._ID)

        return ptnrs

    def print_agent(self):
        #print self.partner_list()
        print "\t%.6d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s" % (self._ID, self._age, self._gender, self._SO, self._DU, self._race, self._HIV_bool,self._incar_bool,self._incar_time, self.partner_list())

    def print_agent_abridge(self):
        print "\t%.6d\t%s\t%s\t%s"%(self._ID, self._gender, self._SO, self._DU)

    def vars(self):
        return "%.6d,%d,%s,%s,%s,%s,%s,%d,%d,%d,%s,%s,%s,%s\n" % (self._ID, self._age, self._gender, self._SO, self._DU, self._race, self._HIV_bool, len(self._partners), self._num_sex_partners, self._timeAlive, self._AIDS_bool, self._tested, self._PrEP_bool, self._incar_bool)

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


class Relationship(object):
    "Class for agent relationships."
    def __init__(self, ID1, ID2, SO, rel_type, duration, initial_agent=False):

        # self._ID is unique ID number used to track each person agent.
        self._ID1 = ID1
        self._ID2 = ID2
        self._ID = self._ID1.get_ID()*100000 + self._ID2.get_ID()
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

    def progress(self, forceKill = False):
        if self._duration <= 0 or forceKill:
            #print self._duration
            agent = self._ID1
            partner = self._ID2

            #print "\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tDeleting relationship between %d and %d" % (agent.get_ID(), partner.get_ID())
            #agent.print_agent()
            #partner.print_agent()
            #self.print_rel()
            agent.unbond(partner, self)
            #agent.unpair(partner)
            #partner.unpair(agent)
            #del self
            return True
            #self._ID1._relationships.remove(self)
            #self._ID2._relationships.remove(self)
        else:
            self._duration = self._duration - 1
            return False

    def delete(self):
        del self

    def get_ID(self):
        return self._IDq

    def print_rel(self):
        #print self.partner_list()
        print "\t%.6d\t%.6d\t%s\t%s\t%d\t%d" % (self._ID1.get_ID(), self._ID2.get_ID(), self._SO, self._rel_type, self._duration, self._total_sex_acts)

    def vars(self):
        return "%.6d,%.6d,%s,%s,%d,%d\n" % (self._ID1.get_ID(), self._ID2.get_ID(), self._SO, self._rel_type, self._duration, self._total_sex_acts)



class Agent_set(Agent):
    """
    Class for agents that contain a "set" of agents from a lower
    hierarchical  level.
    """
    def __init__(self, world, ID, parent = None, numerator = None):
        #Agent.__init__(self, world, ID, initial_agent)

        # _members stores agent set members in a dictionary keyed by ID
        self._ID = ID
        self._members = []
        self._subset = {}

        # _parent_set stores the parent set if this set is a member of an
        # Agent_set class instance. For example, for a set that is a
        # member of a larger set, the _parent_set for that set  would
        # be that larger set.
        self._parent_set = parent
        if numerator:
            self._numerator = numerator
        else:
            self._numerator = self


    def clear_set(self):
        del self._members[:]
        self._subset = {}

    def get_agents(self):
        #return self._members.values()
        return self._members.__iter__()

    def get_agent(self, ID):
        "Returns an agent given the agent's ID"
        return self._members[ID]

    def is_member(self, ID):
        "Returns true if agent is a member of this set"
        return ID in self._members

    def add_agent(self, agent):
        "Adds a new agent to the set."
        #agent.print_agent()
        #if agent in self._members:
        #    raise KeyError("agent %s is already a member of agent set %s"%(agent.get_ID(), self._ID))
        self._members.append(agent)
        #self._members[agent.get_ID()] = agent

        if self._subset: #if subsets exit, check to remove the agent from those lists
            try:
                #self._subset[agent._SO].print_agents()
                self._subset[agent._SO].add_agent(agent)
                #self._members.
            except:
                print "agent %s is already a member of agent set %s"%(agent.get_ID(), self.get_ID())
        # Set the agent's _parent_agent to reflect the parent of this Agent_set
        # instance (self)
        #agent.set_parent_agent(self)

    def remove_agent(self, agent):
        "Removes agent from agent set."

        #ID = agent.get_ID()
        #print "Removing agent %d"%ID

        try:
            self._members.remove(agent)
            #print "agent %s has been removed from agent set %s"%(agent._ID, self.get_ID())
            #self._members.
        except:
            print "agent %s is not a member of agent set %s"%(agent.get_ID(), self.get_ID())
            agent.print_agent()
            #self.print_agents()
            exit(3)

        """Removes agent from SO subset."
        if self._subset: #if subsets exit, check to remove the agent from those lists
            try:
                self._subset[agent._SO]._members.remove(agent)
                #self._members.
            except:
                 print "agent %s is not a member of agent subset %s"%(agent.get_ID(), self._subset[agent._SO].get_ID())
        """

        for tmpS in self.iter_subset():
            #tmpS.print_agents()
            if tmpS.is_member(agent):
                tmpS.remove_agent(agent)
        # Reset the agent's _parent_agent
        #assert agent.get_parent_agent().get_ID() == self.get_ID(), "Removing agent from an Agent_set it does not appear to be assigned to."
        #agent.set_parent_agent(None)
        #print "Removed agent", agent._ID

    def iter_agents(self):
        for agent in self.get_agents():
            yield agent

    def num_members(self):
        return len(self._members)

    def random_agent(self):
        try:
            # tmpA = random.choice(self._members.keys())
            # #print "Returned agent", tmpA._ID
            # return self._members[tmpA]
            return random.choice(self._members)
        except:
            pass
            #print "No agents left", self._members

    def add_subset(self, subset):
        "Adds a new Agent_set to the current sets subset."
        if subset._ID in self._subset:
            raise KeyError("Subset %s is already a member of Agent_set set %s"%(subset.get_ID(), self._ID))
        self._subset[subset._ID] = subset

        # Set the agent's _parent_agent to reflect the parent of this Agent_set
        # instance (self)
        subset.set_parent_agent(self)

    def remove_subset(self, subset):
        "Removes Agent_set to the current sets subset."
        try:
            self._subset.pop(subset.ID)
            #self._members.
        except KeyError:
            raise KeyError("subset %s is not a member of set %s"%(subset.get_ID(), self.get_ID()))

    def iter_subset(self):
        for subset in self._subset.values():
            yield subset

    def set_parent_set(self, master_set):
        self._parent_set = master_set

    def get_parent_set(self):
        return self._parent_set

    def print_agents(self):
        print "\t_____________ %s _____________" % self.get_ID()
        print "\tID\tAge\tGdr\tSO\tDU\tRace\tHIV+\tPtnrs"

        for tmpA in self.iter_agents():
            tmpA.print_agent()
        print "\t______________ END ______________"

    def print_agents_abridged(self):
        print "\t_________ %s _________" % self.get_ID()
        print "\tID\tGdr\tSO\tDU"

        for tmpA in self.iter_agents():
            tmpA.print_agent_abridge()
        print "\t__________ END __________"

    def print_agents_to_file(self, time=None, overWrite="a", filename="Results/Tableau_Agent_Output_File.txt"):
        if overWrite=="a":
            agentList=""
        else:
            agentList="Time,ID,Age,Gdr,SO,DU,Race,HIV+,Ptnrs,LTPtnrs,AliveTime,AIDS,Tested,PrEP,Incar\n"
        for tmpA in self.iter_agents():
            agentList += str(time)+","+tmpA.vars()

        open(filename, overWrite).write(agentList)

    def print_agent_relationshps(self):
        print "\t_____________ %s _____________" % self.get_ID()
        print "\tID1\t\tID2\t\tSO\tTy\tDur\tActs"

        for tmpR in self.iter_agents():
            tmpR.print_rel()
        print "\t______________ END ______________"

    def print_agent_relationships_to_file(self, time=None, overWrite="a"):
        if overWrite=="a":
            agentList=""
        else:
            agentList="Time,ID1,ID2,SO,Ty,Dur,Acts\n"
        for tmpR in self.iter_agents():
            agentList += str(time)+","+tmpR.vars()

        open("Results/Tableau_Rel_Output_File.txt", overWrite).write(agentList)

    def print_subsets(self):
        print "\t__________ %s __________" % self.get_ID()
        print "\tID\t\tMembers\t%"
        for tmpS in self.iter_subset():
            if tmpS.num_members() > 0:
                print "\t%s\t\t%d\t%.2f" % (tmpS._ID, tmpS.num_members(), (1.0*tmpS.num_members()/tmpS._numerator.num_members()))
        print "\t______________ END ______________"

class Agent_Store(object):
    """
    Agent_Store is a class used to store agents who have left for various
    reasons (such as migration) or are in school. It allows triggering their
    return or graduation during a later timestep of the model.
    """
    def __init__(self):
        # self._releases is a dictionary, keyed by timestep, that stores the
        # time at which each agent will be released back to their original
        # parent agent (when they return from school, or from their temporary
        # migration, for example.
        self._releases = {}
        self._parent_dict = {}
        self._stored_agents = []

    def add_agent(self, agent, release_time):
        """
        Adds a new agent to the agent store. Also remove the agent from it's
        parent Agent_set instance.
        """
        if release_time in self._releases:
            self._releases[release_time].append(agent)
        else:
            self._releases[release_time] = [agent]
        self._parent_dict[agent] = agent.get_parent_agent()
        # Store a reference to the agent store with the class instance that is
        # being stored, for easy retrieval later
        agent._store_list.append(self)
        agent.get_parent_agent().remove_agent(agent)
        # Keep a list of the agents stored in this agent_set instance in
        # _stored_agents so that we can easily check whether or now an agent is
        # in an agent store instance, without having to iterate through all the
        # release times.
        self._stored_agents.append(agent)

    def release_agents(self, time):
        # TODO: Make this more general, so it works for households or person
        # agents. Right now it only works for persons since to get the
        # neighborhood ID for tracking we have to call .get_parent_agent twice.
        released_agents = []
        released_agents_dict = {}
        if time in self._releases:
            for agent in self._releases[time]:
                parent_agent = self._parent_dict.pop(agent)
                parent_agent.add_agent(agent)
                agent._store_list.remove(self)
                self._stored_agents.remove(agent)
                neighborhood = parent_agent.get_parent_agent()
                if not neighborhood.get_ID() in released_agents_dict:
                    released_agents_dict[neighborhood.get_ID()] = 0
                released_agents_dict[neighborhood.get_ID()] += 1
                released_agents.append(agent)
            # Remove the now unused releases list for this timestep.
            self._releases.pop(time)
        return released_agents_dict, released_agents

    def in_store(self, agent):
        if agent in self._stored_agents: return True
        else: return False

    def remove_agent(self, agent):
        """
        Remove an agent from the store without releasing it to its original
        location (useful for handling agents who die while away from home).
        """
        self._releases[agent._return_timestep].remove(agent)
        self._parent_dict.pop(agent)
        self._stored_agents.remove(agent)
        agent._store_list.remove(self)

    def __str__(self):
        return 'Agent_Store(%s)'%self._releases