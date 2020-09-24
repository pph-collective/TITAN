from typing import Dict, Optional, ClassVar, Set

from .base_feature import BaseFeature
from ..agent import Agent
from ..population import Population
from ..model import HIVModel
from ..parse_params import ObjMap


class Prep(BaseFeature):

    name = "prep"
    stats = ["prep", "prep_new", "prep_injectable", "prep_oral"]
    """
        PrEP collects the following stats:

        * prep - number of agents with active PrEP
        * prep_new - number of agents who became active PrEP this time step
        * prep_injectable - number of agents on injectable PrEP
        * prep_oral - number of agents on oral PrEP
    """

    # class level attributes to track all Prep agents
    counts: ClassVar[Dict[str, int]] = {}
    new_agents: ClassVar[Set["Agent"]] = set()

    def __init__(self, agent: "Agent"):
        super().__init__(agent)
        # agent level attributes
        self.active = False
        self.adherence = 0
        self.type = ""
        self.load = 0.0
        self.last_dose = 0

    def init_agent(self, pop: "Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        If an agent does not have HIV, is PrEP eligible, and time is at least the prep start time, they are randomly asigned to enroll in PrEP.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        if (
            not self.agent.hiv
            and self.eligible()
            and time >= self.agent.location.params.prep.start_time
            and pop.pop_random.random() < self.agent.location.params.prep.target
        ):
            self.enroll(pop.pop_random)

    def update_agent(self, model: "HIVModel"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `HIVModel.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        If the agent is not hiv and time is at least the prep start time, if the agent is already on PrEP update their PrEP attributes, if the agent isn't on PrEP and is eleigible, initiate PrEP.

        args:
            model: the instance of HIVModel currently being run
        """
        if (
            not self.agent.hiv
            and model.time >= self.agent.location.params.prep.start_time
        ):
            if self.active:
                self.progress(model)
            elif self.eligible():
                self.initiate(model)

    @classmethod
    def add_agent(cls, agent: "Agent"):
        """
        Add an agent to the class (not instance).

        Add agent to the PrEP counts by race, and add the agent to the set of new agents.

        args:
            agent: the agent to add to the class attributes
        """
        # set up if this is the first time being called
        if len(cls.counts) == 0:
            cls.init_class(agent.location.params)

        cls.counts[agent.race] += 1
        cls.new_agents.add(agent)

    @classmethod
    def remove_agent(cls, agent):
        """
        Remove an agent from the class (not instance).

        Decrement the prep counts by race.

        args:
            agent: the agent to remove from the class attributes
        """
        cls.counts[agent.race] -= 1

    @classmethod
    def update_pop(cls, model: "HIVModel"):
        """
        Update the feature for the entire population (class method).

        This is called in `HIVModel.update_all_agents` before agent-level updates are made.

        Resets the tracking set for `new_agents` as this is called before `update_agent`

        args:
            model: the instance of HIVModel currently being run
        """
        cls.new_agents = set()

    def set_stats(self, stats: Dict[str, int]):
        if self.agent in self.new_agents:
            stats["prep_new"] += 1

        if self.active:
            stats["prep"] += 1
            if self.type == "Inj":
                stats["prep_injectable"] += 1
            elif self.type == "Oral":
                stats["prep_oral"] += 1

    # =============== HELPER METHODS ===================

    @classmethod
    def init_class(cls, params: "ObjMap"):
        # internal method to initialize count dictionary
        cls.counts = {race: 0 for race in params.classes.races}

    def initiate(self, model: "HIVModel", force: bool = False):
        """
        Place agents onto PrEP treatment. PrEP treatment assumes that the agent knows their HIV status is negative.

        args:
            model : instance of HIVModel being run
            force : whether to force the agent to enroll instead of using the appropriate algorithm per the prep params
        """
        # Prep only valid for agents not on prep and are HIV negative
        if self.active or self.agent.hiv:
            return

        params = self.agent.location.params

        if len(self.counts) == 0:
            self.init_class(params)

        if force:
            self.enroll(model.run_random)
        else:
            if "Racial" in params.prep.target_model:
                num_prep_agents = self.counts[self.agent.race]
                all_hiv_agents = model.pop.hiv_agents.members
                all_race = {
                    a for a in model.pop.all_agents if a.race == self.agent.race
                }

                num_hiv_agents = len(all_hiv_agents & all_race)
                target_prep = (len(all_race) - num_hiv_agents) * params.demographcis[
                    self.agent.race
                ][self.agent.sex_type].prep.coverage
            else:
                num_prep_agents = sum(self.counts.values())
                target_prep = int(
                    (
                        model.pop.all_agents.num_members()
                        - model.pop.hiv_agents.num_members()
                    )
                    * params.prep.target
                )

            if (
                num_prep_agents < target_prep
                and model.time >= params.prep.start_time
                and self.eligible()
            ):
                self.enroll(model.run_random)

    def enroll(self, rand_gen):
        """
        Enroll an agent in PrEP

        args:
            rand_gen: random number generator
        """
        params = self.agent.location.params

        self.active = True
        self.load = params.prep.peak_load
        self.last_dose = 0

        if (
            rand_gen.random()
            < params.demographics[self.agent.race][self.agent.sex_type].prep.adherence
        ):
            self.adherence = 1
        else:
            self.adherence = 0

        if "Inj" in params.prep.type and "Oral" in params.prep.type:
            if rand_gen.random() < params.prep.lai.prob:
                self.type = "Inj"
            else:
                self.type = "Oral"

        else:
            self.type = params.prep.type[0]

        self.add_agent(self.agent)

    def progress(self, model: "HIVModel", force: bool = False):
        """
        Update agent's PrEP status and discontinue stochastically or if `force` is True

        args:
            model: instance of the HIVModel being run
            force: whether to force discontinuation of PrEP
        """
        if force:
            self.discontinue()  # TO_REVIEW should this just remove the agent from counts, or discontinue? does it depend on type?
            return

        if (
            model.run_random.random()
            < self.agent.location.params.demographics[self.agent.race][
                self.agent.sex_type
            ].prep.discontinue
            and self.type == "Oral"
        ):
            self.discontinue()

        if self.type == "Inj":
            self.update_load()

    def update_load(self):
        """
        Determine and update load of PrEP concentration in agent.
        """
        params = self.agent.location.params

        self.last_dose += 1
        if self.last_dose > params.model.time.steps_per_year:
            self.discontinue()
        else:
            annualized_last_dose = self.last_dose / params.model.time.steps_per_year
            annualized_half_life = params.prep.half_life / 365
            self.load = params.prep.peak_load * (
                (0.5) ** (annualized_last_dose / annualized_half_life)
            )

    def discontinue(self):
        """
        Discontinue PrEP usage
        """
        self.active = False
        self.type = ""
        self.load = 0.0
        self.last_dose = 0

        self.remove_agent(self.agent)

    def eligible(self):
        """
        Determine if an agent is eligible for PrEP

        returns:
            prep
        """
        target_model = self.agent.location.params.prep.target_model
        gender = self.agent.location.params.classes.sex_types[
            self.agent.sex_type
        ].gender

        if (
            self.active
            or self.agent.vaccine.active
            or self.agent.location.params.features.random_trial
        ):
            return False

        all_eligible_models = {"Allcomers", "Racial"}

        if all_eligible_models.intersection(target_model):
            return True

        if "cdc_women" in target_model:
            if gender == "F":
                if self.cdc_eligible():
                    return True

        if "cdc_msm" in target_model:
            if gender == "M" and self.cdc_eligible():
                return True

        if "pwid_sex" in target_model:
            if self.agent.drug_type == "Inj" and self.cdc_eligible():
                return True

        if "pwid" in target_model:
            if self.agent.drug_type == "Inj":
                return True

        if "ssp_sex" in target_model:
            if self.agent.syringe_services.active and self.cdc_eligible():
                return True

        if "ssp" in target_model:
            if self.agent.syringe_services.active:
                return True

        return False

    def cdc_eligible(self) -> bool:
        """
        Determine agent eligibility for PrEP under CDC criteria

        returns:
            cdc eligibility
        """
        if self.agent.is_msm():
            return True

        ongoing_duration = self.agent.location.params.partnership.ongoing_duration
        for rel in self.agent.relationships:
            partner = rel.get_partner(self.agent)
            if rel.duration > ongoing_duration and partner.hiv_dx:
                return True

            if partner.drug_type == "Inj" or partner.is_msm():
                return True

        return False
