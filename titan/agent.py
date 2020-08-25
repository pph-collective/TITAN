#!/usr/bin/env python
# encoding: utf-8

from typing import List, Dict, Set, Optional, Iterable

from .parse_params import ObjMap
from .utils import safe_divide, safe_dist
from .location import Location


class Agent:
    """
    :Purpose:
        This class constructs and represents an agent within the population

    :Input:
        SO : str - sex type flag (HM, HF, MSM)
        age : int - Agents initialization age
        race : str - Race of agent
        DU : str - Drug use flag (IDU, NIDU, NDU)
        location : Location - the location where the agent "lives"
    """

    # class variable for agent creation
    next_agent_id = 0

    @classmethod
    def update_id_counter(cls, last_id):
        cls.next_agent_id = last_id + 1

    def __init__(
        self,
        so: str,
        age: int,
        race: str,
        du: str,
        location: Location,
        id: Optional[int] = None,
    ):
        """
        Initialize an agent based on given properties

        args:
            ID (int) - Unique agent ID
            SO (str) - Sexual orientation flag (HM, HF, MSM)
            age (int) - Agents initialization age
            race (str) - Race of agent
            DU (str) - Drug use flag (Inj, NonInj, None)

        returns:
            None
        """
        # self.id is unique ID number used to track each person agent.
        if id is not None:
            self.id = id
        else:
            self.id = self.next_agent_id

        self.update_id_counter(self.id)

        # agent properties
        self.sex_type = so
        self.age = age
        self.age_bin = 0
        self.race = race
        self.drug_type = du
        self.location = location

        if self.drug_type == "Inj":
            self.population = "PWID"
        else:
            self.population = self.sex_type

        self.msmw = False
        self.sex_role = "versatile"

        # agent-partner params
        self.relationships: Set[Relationship] = set()
        self.partners: Dict[str, Set] = {}
        self.mean_num_partners: Dict[str, int] = {}
        self.target_partners: Dict[str, int] = {}

        # agent STI params
        self.hiv = False
        self.hiv_time = 0
        self.hiv_dx = False
        self.aids = False

        # agent treatment params
        self.haart = False
        self.haart_time = 0
        self.haart_adherence = 0
        self.ssp = False
        self.prep = False
        self.prep_adherence = 0
        self.prep_reason: List[str] = []
        self.intervention_ever = False
        self.random_trial_enrolled = False
        self.vaccine = False
        self.vaccine_time = 0
        self.vaccine_type = ""
        self.partner_traced = False
        self.trace_time = 0
        self.prep_awareness = False
        self.prep_opinion = 0.0
        self.prep_type = ""
        self.pca = False
        self.pca_suitable = False

        # PrEP pharmacokinetics
        self.prep_load = 0.0
        self.prep_last_dose = 0

        # agent high risk params
        self.high_risk = False
        self.high_risk_time = 0
        self.high_risk_ever = False

        # agent incarcartion params
        self.incar = False
        self.incar_time = 0

    def __str__(self):
        """
        String formatting of agent object

        returns:
            String formatted tab-deliminated agent properties
        """
        return (
            f"\t{self.id}\t{self.age}\t{self.sex_type}\t{self.drug_type}\t"
            f"{self.race}\t{self.hiv}"
        )

    def __repr__(self):
        """
        Repr formatting of agent object

        returns:
            ID (str) - agent ID as str
        """
        return str(self.id)

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id

    def iter_partners(self) -> Iterable["Agent"]:
        for partner_set in self.partners.values():
            for partner in partner_set:
                yield partner

    def has_partners(self) -> bool:
        return any(self.iter_partners())

    def get_acute_status(self, acute_time_period) -> bool:
        """
        :Purpose:
            Get acute status of agent at time period
        :Input:
            None
        :Output:
            acute_status : bool
        """
        hiv_t = self.hiv_time

        if acute_time_period >= hiv_t > 0:
            return True
        else:
            return False

    def cdc_eligible(self, ongoing_duration: int, sex_def) -> bool:
        """
        Determine agent eligibility for PrEP under CDC criteria

        args:
            ongoing_duration: Number of time steps above which a relationship is considered ongoing [params.partnership.ongoing_duration]
            sex_def: definition of agent's sex type [param.classes.sex_type]

        returns:
            cdc eligibility
        """
        eligible = False
        for rel in self.relationships:
            if rel.duration > ongoing_duration:
                partner = rel.get_partner(self)
                if partner.drug_type == "Inj":
                    eligible = True
                    self.prep_reason.append("PWID")
                if partner.hiv_dx:
                    eligible = True
                    self.prep_reason.append("HIV test")
                if partner.msmw:
                    eligible = True
                    self.prep_reason.append("MSMW")

        if sex_def.gender == "M":
            if "MSM" in sex_def.sleeps_with:
                eligible = True

        return eligible

    def prep_eligible(self, target_model: set, ongoing_duration: int, sex_def) -> bool:
        """
            Determine if an agent is eligible for PrEP under target model

        args:
            target_model: Model of prep allocation [params.prep.target_model]
            ongoing_duration: Number of time steps above which a relationship is considered ongoing [params.partnership.ongoing_duration]
            sex_def: definition of agent's sex type [param.classes.sex_type]

        returns
            eligibility
        """
        # if agent is already on prep, not eligible to enroll
        if self.prep or self.vaccine:
            return False

        target_model = set(target_model)

        if target_model.intersection({"Allcomers", "Racial"}):
            return True

        if "cdc_women" in target_model:
            if self.sex_type == "HF":
                if self.cdc_eligible(ongoing_duration, sex_def):
                    return True

        if "cdc_msm" in target_model:
            if self.sex_type in ("MSM", "MTF") and self.cdc_eligible(
                ongoing_duration, sex_def
            ):
                return True

        if "pwid_sex" in target_model:
            if self.drug_type == "Inj" and self.cdc_eligible(ongoing_duration, sex_def):
                return True

        if "pwid" in target_model:
            if self.drug_type == "Inj":
                return True

        if "ssp_sex" in target_model:
            if self.ssp and self.cdc_eligible(ongoing_duration, sex_def):
                return True

        if "ssp" in target_model:
            if self.ssp:
                return True

        return False

    def enroll_prep(self, rand_gen):
        params = self.location.params
        self.prep = True
        self.prep_load = params.prep.peak_load
        self.prep_last_dose = 0

        if (
            rand_gen.random()
            < params.demographics[self.race][self.sex_type].prep.adherence
        ):
            self.prep_adherence = 1
        else:
            self.prep_adherence = 0

        # set PrEP load and dosestep for PCK
        if "Inj" in params.prep.type and "Oral" in params.prep.type:
            if rand_gen.random() < params.prep.lai.prob:
                self.prep_type = "Inj"
            else:
                self.prep_type = "Oral"
        else:
            self.prep_type = params.prep.type[0]

    def update_prep_load(self):
        """
        :Purpose:
            Determine and update load of PrEP concentration in agent.

        :Input:
            none

        :Output:
            none
        """
        params = self.location.params
        # N(t) = N0 (0.5)^(t/t_half)
        self.prep_last_dose += 1
        if self.prep_last_dose > params.model.time.steps_per_year:
            self.prep_load = 0.0
            self.prep = False
            self.prep_reason = []
            self.prep_type = ""
            self.prep_last_dose = 0
        else:
            annualized_last_dose = (
                self.prep_last_dose / params.model.time.steps_per_year
            )
            annualized_half_life = params.prep.half_life / 365
            self.prep_load = params.prep.peak_load * (
                (0.5) ** (annualized_last_dose / annualized_half_life)
            )

    def vaccinate(self, vax: str):
        """
        :Purpose:
            Vaccinate an agent and update relevant fields.

        :Input:
            vax - str : Vaccine type

        :Output:
            none
        """
        self.vaccine = True
        self.vaccine_type = vax
        self.vaccine_time = 1

    def get_partners(self, bond_types: Optional[Iterable[str]] = None) -> Set["Agent"]:
        if bond_types:
            partners = set()
            for bond in bond_types:
                partners.update(self.partners[bond])
        else:
            partners = {partner for partner in self.iter_partners()}

        return partners

    def get_num_partners(self, bond_types: Optional[Iterable[str]] = None) -> int:
        return len(self.get_partners(bond_types))

    def get_number_of_sex_acts(self, rand_gen) -> int:
        """
        :Purpose:
            Number of sexActs an agent has done.

        :Input:
            rand_gen : random number generator (e.g. self.runRandom in model)

        :Output:
            number_sex_act : int
        """
        freq_params = self.location.params.partnership.sex.frequency
        if freq_params.type == "bins":
            rv = rand_gen.random()

            for i in range(1, 6):
                p = freq_params.bins[i].prob
                if rv <= p:
                    break

            min_frequency = freq_params.bins[i].min
            max_frequency = freq_params.bins[i].max
            return rand_gen.randint(min_frequency, max_frequency)

        elif freq_params.type == "distribution":
            return round(safe_dist(freq_params.distribution, rand_gen))
        else:
            raise Exception("Sex acts must be defined as bin or distribution")


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
        :Purpose:
            Constructor for a Relationship

        :Input:
            :ID1: first agent
            :ID2: second agent
            :SO: Orientation of relationship
            :bond_type: (future feature) - sex bond or idu bond or both
            :duration: length of relationship
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
        # TODO MAKE THIS INCREMENT WITH passed IDs

        # Relationship properties
        self.duration = duration
        self.total_duration = duration
        self.total_sex_acts = 0
        self.bond_type = bond_type

        self.bond()

    def __eq__(self, other):
        return self.id == other.id

    def __hash__(self):
        return self.id

    def progress(self, force: bool = False):
        if self.duration <= 0 or force:
            self.unbond()
            return True
        else:
            self.duration -= 1
            return False

    def bond(self):
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
        self.agent1.relationships.add(self)
        self.agent2.relationships.add(self)

        # Pair agent with partner and partner with agent
        self.agent1.partners[self.bond_type].add(self.agent2)
        self.agent2.partners[self.bond_type].add(self.agent1)

    def unbond(self):
        """
        Unbond two agents. Removes relationship from relationship lists.
        Removes partners in each others' partner list.

        args:
            agent: Agent - former partner of partner
            partner: Agent - former partner of agent
        returns:
            None
        """

        # Remove relationship to relationships list for each agent
        self.agent1.relationships.remove(self)
        self.agent2.relationships.remove(self)

        # Unpair agent with partner and partner with agent
        self.agent1.partners[self.bond_type].remove(self.agent2)
        self.agent2.partners[self.bond_type].remove(self.agent1)

    def get_partner(self, agent: "Agent") -> "Agent":
        if agent == self.agent1:
            return self.agent2
        else:
            return self.agent1

    def __str__(self):
        return (
            f"\t{self.id}\t{self.agent1.id}\t{self.agent2.id}\t{self.duration}\t"
            f"{self.bond_type} "
        )

    def __repr__(self):
        return str(self.id)


class AgentSet:
    """
    Class for agents that contain a "set" of agents from a lower
    hierarchical  level.
    """

    def __init__(self, id: str, parent: "AgentSet" = None):
        # _members stores agent set members in a dictionary keyed by ID
        self.id = id
        self.members: Set[Agent] = set()
        self.subset: Dict[str, AgentSet] = {}

        # _parent_set stores the parent set if this set is a member of an
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
        self.members: Set[Agent] = set()
        self.subset: Dict[str, str] = {}

    def __iter__(self):
        return self.members.__iter__()

    def __contains__(self, item):
        return self.members.__contains__(item)

    def is_member(self, agent: Agent):
        """Returns true if agent is a member of this set"""
        return agent in self.members

    # adding trickles up
    def add_agent(self, agent: Agent):
        """Adds a new agent to the set."""
        self.members.add(agent)

        if self.parent_set is not None:
            self.parent_set.add_agent(agent)

    # removing trickles down
    def remove_agent(self, agent: Agent):
        """Removes agent from agent set."""
        if agent in self.members:
            self.members.remove(agent)

        for subset in self.iter_subset():
            subset.remove_agent(agent)

    def num_members(self) -> int:
        return len(self.members)

    def add_subset(self, subset: "AgentSet"):
        """Adds a new AgentSet to the current sets subset."""
        if subset.id not in self.subset:
            self.subset[subset.id] = subset

    def iter_subset(self):
        for subset in list(self.subset.values()):
            yield subset

    def print_subsets(self):
        print(f"\t__________ {self.id} __________")
        print("\tID\t\tN\t\t%")
        for set in self.iter_subset():
            print(
                "\t{:6}\t\t{:5}\t\t{:.2}".format(
                    set.id,
                    set.num_members(),
                    safe_divide(set.num_members(), set.parent_set.num_members()),
                )
            )
            for subset in set.iter_subset():
                print(
                    "\t{:4}\t\t{:5}\t\t{:.2}".format(
                        subset.id,
                        subset.num_members(),
                        safe_divide(
                            subset.num_members(), subset.parent_set.num_members()
                        ),
                    )
                )
        print("\t______________ END ______________")
