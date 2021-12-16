#!/usr/bin/env python
# encoding: utf-8

from typing import Dict, Set, Optional, Iterator, Iterable

from .utils import (
    safe_divide,
    safe_dist,
    get_independent_bin,
    safe_random_choice,
    safe_random_int,
)
from .location import Location
from . import features
from . import exposures


class Agent:
    """
    This class constructs and represents an agent within the population
    """

    # class variable for agent creation
    next_agent_id = 0

    @classmethod
    def update_id_counter(cls, last_id):
        cls.next_agent_id = last_id + 1

    def __init__(
        self,
        sex_type: str,
        age: int,
        race: str,
        drug_use: str,
        location: Location,
        id: Optional[int] = None,
    ) -> None:
        """
        Initialize an agent based on given properties

        Args:
            id: Unique agent ID
            sex_type: Name of defined sex type (e.g. MSM) [params.classes.sex_types]
            age: Agents initialization age
            race: Race of agent [params.classes.races]
            drug_use: Drug use flag [params.classes.drug_types]
        """
        # self.id is unique ID number used to track each person agent.
        if id is not None:
            self.id = id
        else:
            self.id = self.next_agent_id

        self.update_id_counter(self.id)

        # agent properties
        self.sex_type = sex_type
        self.age = age
        self.race = race
        self.drug_type = drug_use
        self.location = location
        self.component = "-1"  # updated after relationships created

        self.sex_role = "versatile"

        # agent-partner params
        self.relationships: Set[Relationship] = set()
        self.partners: Dict[str, Set] = {}
        self.mean_num_partners: Dict[str, int] = {}
        self.target_partners: Dict[str, int] = {}

        # agent exposures params
        # model features
        for exposure in exposures.BaseExposure.__subclasses__():
            setattr(self, exposure.name, exposure(self))

        # model features
        for feature in features.BaseFeature.__subclasses__():
            setattr(self, feature.name, feature(self))

    def __str__(self) -> str:
        """
        String formatting of agent object

        returns:
            String formatted tab-deliminated agent properties
        """
        return (
            f"\t{self.id}\t{self.age}\t{self.sex_type}\t{self.drug_type}\t"  # type: ignore[attr-defined]
            f"{self.race}\t{self.hiv.active}"  # type: ignore[attr-defined]
        )

    def __repr__(self) -> str:
        """
        Repr formatting of agent object

        returns:
            agent ID as str
        """
        return str(self.id)

    def __eq__(self, other) -> bool:
        return self.id == other.id

    def __ne__(self, other) -> bool:
        return self.id != other.id

    def __hash__(self) -> int:
        return self.id

    def iter_partners(self) -> Iterator["Agent"]:
        """
        Get an iterator over an agent's partners

        returns:
            iterator of agent partners
        """
        for partner_set in self.partners.values():
            for partner in partner_set:
                yield partner

    def has_partners(self) -> bool:
        """
        Determine whether an agent has any partners

        returns:
            whether an agent has at least one partner
        """
        return any(self.iter_partners())

    def is_msm(self) -> bool:
        """
        Determine whether an agent is a man who can have sex with men

        returns:
            if agent is MSM
        """
        sex_dict = self.location.params.classes.sex_types
        if sex_dict[self.sex_type].gender != "M":
            return False

        for sex_type in sex_dict[self.sex_type].sleeps_with:
            if sex_dict[sex_type].gender == "M":
                return True
        return False

    def get_partners(self, bond_types: Optional[Iterable[str]] = None) -> Set["Agent"]:
        """
        Get all of an agents partners or those with specific bond types

        args:
            bond_types: list of bond types which will filter the partners, otherwise all partners returned

        returns:
            set of agent's partners
        """
        if bond_types:
            partners = set()
            for bond in bond_types:
                partners.update(self.partners[bond])
        else:
            partners = {partner for partner in self.iter_partners()}

        return partners

    def get_num_partners(self, bond_types: Optional[Iterable[str]] = None) -> int:
        """
        Get the number of partners an agent has, optionally filtered by bond type

        args:
            bond_types: list of bond types which will filter the partners, otherwise total number of partners returned

        returns:
            the number of partners the agent has
        """
        return len(self.get_partners(bond_types))


class Relationship:
    """Class for agent relationships."""

    # class variable for relationship creation
    next_rel_id = 0

    @classmethod
    def update_id_counter(cls, last_id):
        cls.next_rel_id = last_id + 1

    def __init__(
        self,
        agent1: Agent,
        agent2: Agent,
        duration: int,
        bond_type: str,
        id: Optional[int] = None,
    ):
        """
        Constructor for a Relationship

        args:
            agent1: first agent
            agent2: second agent
            duration: target duration of relationship
            bond_type: type of bond for the relationship [params.classes.bond_types]
            id: unique identifier
        """
        # make sure these agents can be in a relationship
        assert agent1 != agent2, "Cannot create relationship with same agent"
        for rel in agent1.relationships:
            assert agent2 != rel.get_partner(agent1), "Agents already partnered!"

        # self.id is unique ID number used to track each person agent.
        self.agent1 = agent1
        self.agent2 = agent2

        if id is not None:
            self.id = id
        else:
            self.id = self.next_rel_id

        self.update_id_counter(self.id)

        # Relationship properties
        self.duration = duration
        self.total_duration = duration
        self.bond_type = bond_type

        self.bond()

    def __eq__(self, other) -> bool:
        return self.id == other.id

    def __ne__(self, other) -> bool:
        return self.id != other.id

    def __hash__(self) -> int:
        return self.id

    def progress(self, force: bool = False) -> bool:
        """
        Progress a relationship to the next time step (decrementing remaining target duration), or end a relationship if the duration is 0 or if `force` is set to `True`

        args:
            force: whether to force the relationship to end
        """
        if self.duration <= 0 or force:
            self.unbond()
            return True
        else:
            self.duration -= 1
            return False

    def bond(self) -> None:
        """
        Bond two agents. Adds the relationship to each agent's `relationships` set, then adds each partner to the others' partner list.
        """

        # Append relationship to relationships list for each agent
        self.agent1.relationships.add(self)
        self.agent2.relationships.add(self)

        # Pair agent with partner and partner with agent
        self.agent1.partners[self.bond_type].add(self.agent2)
        self.agent2.partners[self.bond_type].add(self.agent1)

    def unbond(self):
        """
        Unbond two agents. Removes relationship from relationship sets.
        Removes partners in each others' partner list.
        """

        # Remove relationship to relationships list for each agent
        self.agent1.relationships.remove(self)
        self.agent2.relationships.remove(self)

        # Unpair agent with partner and partner with agent
        self.agent1.partners[self.bond_type].remove(self.agent2)
        self.agent2.partners[self.bond_type].remove(self.agent1)

    def get_partner(self, agent: "Agent") -> "Agent":
        """
        Given an agent in the relationship, return the other agent

        args:
            agent: one of the agents in the relationship

        returns:
            the agent's partner
        """
        if agent == self.agent1:
            return self.agent2
        elif agent == self.agent2:
            return self.agent1
        else:
            raise ValueError("Agent must be in this relationship")

    def get_number_of_sex_acts(self, rand_gen) -> int:
        """
        Number of sex acts in the relationship during the time step.

        args:
            rand_gen: np random number generator (e.g. self.run_random in model)

        returns:
            number of sex acts
        """
        agent = safe_random_choice([self.agent1, self.agent2], rand_gen)
        freq_params = agent.location.params.partnership.sex.frequency[self.bond_type]

        if freq_params.type == "bins":
            i = get_independent_bin(rand_gen, freq_params.bins)
            return safe_random_int(
                freq_params.bins[i].min, freq_params.bins[i].max, rand_gen
            )

        elif freq_params.type == "distribution":
            return round(safe_dist(freq_params.distribution, rand_gen))

        else:
            raise Exception("Sex acts must be defined as bin or distribution")

    def __str__(self):
        return (
            f"\t{self.id}\t{self.agent1.id}\t{self.agent2.id}\t{self.duration}\t"
            f"{self.bond_type} "
        )

    def __repr__(self):
        return str(self.id)


class AgentSet:
    """
    Container for agents into heirarchical sets (e.g. all_agents > hiv_agents)
    """

    def __init__(
        self,
        id: str,
        parent: Optional["AgentSet"] = None,
    ):
        """
        Constructor of an AgentSet

        args:
            id: name of the set
            parent: the set this set is a subset of
        """
        # members stores agent set members in a dictionary keyed by ID
        self.id = id
        self.members: Set[Agent] = set()
        self.subset: Dict[str, AgentSet] = {}

        # parent_set stores the parent set if this set is a member of an
        # AgentSet class instance. For example, for a set that is a
        # member of a larger set, the _parent_set for that set  would
        # be that larger set.
        self.parent_set = parent
        if parent:
            parent.add_subset(self)

    def __repr__(self):
        return self.id

    def __str__(self):
        return self.id

    def clear_set(self):
        """
        Clears a set of any members and subsets
        """
        self.members: Set[Agent] = set()
        self.subset: Dict[str, str] = {}

    def __iter__(self) -> Iterator[Agent]:
        """
        Iterate over the memers in the set

        example:
            ```py
            for agent in agent_set:
                # do something with agent
            ```

        returns:
            iterator over member agents
        """
        return self.members.__iter__()

    def __contains__(self, item) -> bool:
        """
        Is an agent a member of this agent set

        example:
            ```py
            if agent in agent_set:
                # do something
            ```

        returns:
            whether agent is part of set
        """
        return self.members.__contains__(item)

    # adding trickles up
    def add_agent(self, agent: Agent) -> None:
        """
        Adds an agent to the set and any parent sets

        args:
            agent: agent to add
        """
        self.members.add(agent)

        if self.parent_set is not None:
            self.parent_set.add_agent(agent)

    # removing trickles down
    def remove_agent(self, agent: Agent) -> None:
        """
        Removes agent from agent set if they are a member of the set.  Also removes the agent from any subsets.

        args:
            agent: agent to remove
        """
        if agent in self.members:
            self.members.remove(agent)

        for subset in self.iter_subset():
            subset.remove_agent(agent)

    def num_members(self) -> int:
        """
        Number of members in the set

        returns:
            number of members
        """
        return len(self.members)

    def add_subset(self, subset: "AgentSet") -> None:
        """
        Adds a new AgentSet to the current sets subset.

        args:
            subset: subset to add to this set
        """
        if subset.id not in self.subset:
            self.subset[subset.id] = subset

    def iter_subset(self) -> Iterator["AgentSet"]:
        """
        Iterate over the subsets of this agent set

        returns:
            iterator of agent sets
        """
        for subset in list(self.subset.values()):
            yield subset

    def print_subsets(self, printer=print):
        """
        Pretty print the subsets of this agent set
        """
        lines = []

        lines.append(f"\t__________ {self.id} __________")
        # lines.append("\tID\t\tN\t\t%")
        lines.append("\t{:^6}\t\t{:^5}\t\t{:^4}".format("ID", "N", "%"))
        for set in self.iter_subset():
            lines.append(
                "\t{:^6}\t\t{:^5}\t\t{:.2}".format(
                    set.id,
                    set.num_members(),
                    safe_divide(set.num_members(), set.parent_set.num_members()),
                )
            )
            for subset in set.iter_subset():
                lines.append(
                    "\t{:4}\t\t{:5}\t\t{:.2}".format(
                        subset.id,
                        subset.num_members(),
                        safe_divide(
                            subset.num_members(), subset.parent_set.num_members()
                        ),
                    )
                )
        lines.append("\t______________ END ______________")
        printer("\n".join(lines))
