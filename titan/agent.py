#!/usr/bin/env python
# encoding: utf-8

from typing import Sequence, List, Dict, Optional

from dotmap import DotMap  # type: ignore


class Agent:
    "Class for agent objects."

    # class variable for agent creation
    next_agent_id = 0

    @classmethod
    def update_id_counter(cls):
        cls.next_agent_id += 1

    def __init__(self, SO: str, age: int, race: str, DU: str):
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
        self.id = self.next_agent_id
        self.update_id_counter()

        # agent properties
        self.so = SO  # REVIEWED split this out into gender and sleeps_with
        self.age = age
        self.age_bin = 0
        self.race = race
        self.drug_use = DU

        self.msmw = False

        # agent-partner params
        self.relationships: List[Relationship] = []
        self.partners: List[Agent] = []
        self.mean_num_partners = 0

        # agent STI params
        self.hiv = False
        self.hiv_time = 0
        self.hiv_dx = False
        self.aids = False

        # agent treatment params
        self.haart = False
        self.haart_time = 0
        self.haart_adherence = 0
        self.sne = False
        self.prep = False
        self.prep_adherence = 0
        self.prep_reason: List[str] = []
        self.intervention_ever = False
        self.vaccine = False
        self.vaccine_time = 0
        self.vaccine_type = ""
        self.partner_traced = False
        self.trace_time = 0
        self.awareness = False
        self.opinion = 0.0
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
        self.incar_ever = False  # TO_REVIEW not really used
        self.incar_time = 0

    def __str__(self):
        """
        String formatting of agent object

        returns:
            String formatted tab-deliminated agent properties
        """
        return "\t%.6d\t%d\t%s\t%s\t%s\t%s" % (
            self.id,
            self.age,
            self.so,
            self.drug_use,
            self.race,
            self.hiv,
        )

    def __repr__(self):
        """
        Repr formatting of agent object

        returns:
            ID (str) - agent ID as str
        """
        return str(self.id)

    def __eq__(self, other):
        if not isinstance(other, Agent):
            return NotImplemented
        return self.id == other.id

    def __hash__(self):
        return hash(self.id)

    def partner_list(self):
        """
        Return the list of partners for this agent

        returns:
            _partners (list) - list of partners
        """
        ptnrs = list()
        if self.partners is not None:
            for partner in self.partners:
                ptnrs.append(partner.id)

        return ptnrs

    def get_acute_status(self) -> bool:
        """
        :Purpose:
            Get acute status of agent at time period
        :Input:
            None
        :Output:
            acute_status : bool
        """
        acute_time_period = 2
        hiv_t = self.hiv_time

        if hiv_t <= acute_time_period and hiv_t > 0:
            return True
        else:
            return False

    def PrEP_eligible(self, target_model: str) -> bool:
        """
        :Purpose:
            Determine if an agent is eligible for PrEP

        :Input:
            time : int

        :Output:
            eligibility : bool
        """
        # if agent is already on prep, not eligible to enroll
        if self.prep:
            return False

        eligible = False
        if target_model in ("Allcomers", "Racial"):
            eligible = True
        elif target_model == "CDCwomen":
            if self.so == "HF":
                for rel in set(self.relationships):
                    partner = rel.get_partner(self)
                    if rel.duration > 1:
                        if partner.drug_use == "Inj":
                            eligible = True
                            self.prep_reason.append("PWID")
                        if partner.hiv_dx:
                            eligible = True
                            self.prep_reason.append("HIV test")
                        if partner.msmw:
                            eligible = True
                            self.prep_reason.append("MSMW")
        elif target_model == "CDCmsm":
            if self.so == "MSM":
                for rel in self.relationships:
                    partner = rel.get_partner(self)

                    if rel.duration > 1:
                        if partner.hiv_dx or self.mean_num_partners > 1:
                            eligible = True
        elif target_model == "MSM":
            if self.so in ("MSM", "MTF"):
                eligible = True

        return eligible

    def update_prep_load(self, params: DotMap):
        """
        :Purpose:
            Determine and update load of PrEP concentration in agent.

        :Input:
            none

        :Output:
            none
        """
        # N(t) = N0 (0.5)^(t/t_half)
        self.prep_last_dose += 1
        if self.prep_last_dose > 12:
            self.prep_load = 0.0
        else:
            self.prep_load = params.prep.peak_load * (
                (0.5) ** (self.prep_last_dose / (params.prep.half_life / 30))
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

    def get_transmission_probability(self, interaction: str, params: DotMap) -> float:
        """ Decriptor
        :Purpose:
            Determines the probability of a transmission event based on interaction type.

        :Input:
            interaction : str - "NEEDLE" or "SEX"

        :Output:
            probability : float
        """
        # Logic for if needle or sex type interaction
        p: float
        if interaction == "NEEDLE":
            p = params.partnership.needle.transmission[self.haart_adherence].prob
        elif interaction == "SEX":
            p = params.partnership.sex.transmission[self.so][self.haart_adherence].prob

        # Scaling parameter for acute HIV infections
        if self.get_acute_status():
            p *= params.calibration.acute

        # Scaling parameter for positively identified HIV agents
        if self.hiv_dx:
            p *= 1 - params.calibration.risk_reduction.transmission

        # Tuning parameter for ART efficiency
        if self.haart:
            p *= params.calibration.risk_reduction.haart

        # Racial calibration parameter to attain proper race incidence disparity
        if self.race == "BLACK":
            p *= params.calibration.race_transmission

        # Scaling parameter for per act transmission.
        p *= params.calibration.transmission

        return p

    def get_number_of_sex_acts(self, rand_gen, params: DotMap) -> int:
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

        for i in range(1, 6):
            p = params.partnership.sex.frequency[i].prob
            if rv <= p:
                min = params.partnership.sex.frequency[i].min
                max = params.partnership.sex.frequency[i].max
                return rand_gen.randrange(min, max, 1)

        # fallthrough is last i
        min = params.partnership.sex.frequency[i].min
        max = params.partnership.sex.frequency[i].max
        return rand_gen.randrange(min, max, 1)


class Relationship:
    "Class for agent relationships."

    # class variable for relationship creation
    next_rel_id = 0

    @classmethod
    def update_id_counter(cls):
        cls.next_rel_id += 1

    def __init__(self, agent1: Agent, agent2: Agent, duration: int, bond_type: str):
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
        assert (
            agent1 not in agent2.partners and agent2 not in agent1.partners
        ), "Agent's already partnered!"

        # self.id is unique ID number used to track each person agent.
        self.agent1 = agent1
        self.agent2 = agent2
        self.id = self.next_rel_id
        self.update_id_counter()

        # Relationship properties
        self.duration = duration
        self.total_duration = duration
        self.total_sex_acts = 0
        self.bond_type = bond_type

        self.bond(agent1, agent2)

    def __eq__(self, other):
        if not isinstance(other, Relationship):
            return NotImplemented
        return self.id == other.id

    def __hash__(self):
        return hash(self.id)

    def progress(self, force: bool = False):
        if self.duration <= 0 or force:
            self.unbond(self.agent1, self.agent2)
            return True
        else:
            self.duration = self.duration - 1
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
        agent.relationships.append(self)
        partner.relationships.append(self)

        # Pair agent with partner and partner with agent
        agent.partners.append(partner)
        partner.partners.append(agent)

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
        agent.relationships.remove(self)
        partner.relationships.remove(self)

        # Unpair agent with partner and partner with agent
        agent.partners.remove(partner)
        partner.partners.remove(agent)

    def get_partner(self, agent: "Agent") -> "Agent":
        if agent == self.agent1:
            return self.agent2
        else:
            return self.agent1

    def __repr__(self):
        return "\t%.6d\t%.6d\t%s\t%s\t%d\t%d" % (
            self.agent1.id,
            self.agent2.id,
            self.duration,
            self.total_sex_acts,
        )


class AgentSet:
    """
    Class for agents that contain a "set" of agents from a lower
    hierarchical  level.
    """

    def __init__(
        self, id: str, parent: "AgentSet" = None, numerator: "AgentSet" = None
    ):
        # _members stores agent set members in a dictionary keyed by ID
        self.id = id
        self.members: List[Agent] = []
        self.subset: Dict[str, AgentSet] = {}

        # _parent_set stores the parent set if this set is a member of an
        # AgentSet class instance. For example, for a set that is a
        # member of a larger set, the _parent_set for that set  would
        # be that larger set.
        self.parent_set = parent
        if parent:
            parent.add_subset(self)
        if numerator:
            self.numerator = numerator
        else:
            self.numerator = self

    def __repr__(self):
        return self.id

    def __str__(self):
        return self.id

    def clear_set(self):
        self.members: List[Agent] = []
        self.subset: Dict[str, str] = {}

    def __iter__(self):
        return self.members.__iter__()

    def is_member(self, agent: Agent):
        "Returns true if agent is a member of this set"
        return agent in self.members

    # adding trickles up
    def add_agent(self, agent: Agent):
        "Adds a new agent to the set."
        if not self.is_member(agent):
            self.members.append(agent)

            if self.parent_set is not None:
                self.parent_set.add_agent(agent)

    # removing trickles down
    def remove_agent(self, agent: Agent):
        "Removes agent from agent set."
        if self.is_member(agent):
            self.members.remove(agent)

        for tmpS in self.iter_subset():
            tmpS.remove_agent(agent)

    def num_members(self) -> int:
        return len(self.members)

    def add_subset(self, subset: "AgentSet"):
        "Adds a new AgentSet to the current sets subset."
        if subset.id not in self.subset:
            self.subset[subset.id] = subset

    def iter_subset(self):
        for subset in list(self.subset.values()):
            yield subset

    def print_subsets(self):
        print("\t__________ %s __________" % self.id)
        print("\tID\t\tN\t\t%")
        for tmpS in self.iter_subset():
            if tmpS.num_members() > 0:
                print(
                    "\t%-6s\t%-5d\t%.2f"
                    % (
                        tmpS.id,
                        tmpS.num_members(),
                        (1.0 * tmpS.num_members() / tmpS.numerator.num_members()),
                    )
                )
            for tmpSS in tmpS.iter_subset():
                if tmpSS.num_members() > 0:
                    print(
                        "\t-%-4s\t%-5d\t%.2f"
                        % (
                            tmpSS.id,
                            tmpSS.num_members(),
                            (1.0 * tmpSS.num_members() / tmpSS.numerator.num_members()),
                        )
                    )
        print("\t______________ END ______________")
