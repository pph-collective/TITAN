from typing import Dict, ClassVar, Optional

import numpy as np  # type: ignore

from . import base_feature
from .. import agent
from .. import population
from .. import model
from ..parse_params import ObjMap


class Prep(base_feature.BaseFeature):

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

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)
        # agent level attributes
        self.active = False
        self.adherent = False
        self.type = ""
        self.time = None
        self.last_dose_time: Optional[int] = None

    @classmethod
    def init_class(cls, params: "ObjMap"):
        """
        Initialize the counts dictionary for the races in the model.

        args:
            params: the population params
        """
        cls.counts = {race: 0 for race in params.classes.races}

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        If an agent does not have HIV, is PrEP eligible, and time is at least the prep start time, they are randomly asigned to enroll in PrEP.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        params = self.agent.location.params
        if not self.agent.hiv and self.eligible() and time >= params.prep.start_time:
            if params.prep.target_as_prob:
                if "Racial" in params.prep.target_model:
                    if (
                        pop.pop_random.random()
                        < params.demographics[self.agent.race][
                            self.agent.sex_type
                        ].prep.init
                    ):
                        self.enroll(pop.pop_random, time)
                elif pop.pop_random.random() < params.prep.init:
                    self.enroll(pop.pop_random, time)
            elif "Racial" in params.prep.target_model:
                if (
                    pop.pop_random.random()
                    < params.demographics[agent.race][agent.sex_type].prep.target
                ):
                    self.enroll(pop.pop_random, time)
            elif pop.pop_random.random() < params.prep.target:
                self.enroll(pop.pop_random, time)

    def update_agent(self, model: "model.HIVModel"):
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
    def add_agent(cls, agent: "agent.Agent"):
        """
        Add an agent to the class (not instance).

        Add agent to the PrEP counts by race, and add the agent to the set of new agents.

        args:
            agent: the agent to add to the class attributes
        """
        # set up if this is the first time being called
        cls.counts[agent.race] += 1

    @classmethod
    def remove_agent(cls, agent):
        """
        Remove an agent from the class (not instance).

        Decrement the prep counts by race.

        args:
            agent: the agent to remove from the class attributes
        """
        cls.counts[agent.race] -= 1

    def set_stats(self, stats: Dict[str, int], time: int):
        if self.active:
            stats["prep"] += 1

            if self.time == time:
                stats["prep_new"] += 1

            if self.type == "Inj":
                stats["prep_injectable"] += 1
            elif self.type == "Oral":
                stats["prep_oral"] += 1

    def get_acquisition_risk_multiplier(self, time: int, interaction_type: str):
        """
        Get a multiplier for how prep reduces risk of HIV acquisition.

        By default, returns 1.0

        args:
            time: the current model time step
            interaction_type: The type of interaction where the agent could acquire HIV (e.g. 'sex', 'injection' - from [params.classes.interaction_types])
        """
        if self.active and self.last_dose_time is not None:
            params = self.agent.location.params
            if self.type == "Oral":
                if self.adherent:
                    return 1.0 - params.prep.efficacy.adherent
                else:
                    return 1.0 - params.prep.efficacy.non_adherent
            elif self.type == "Inj" and not self.adherent:
                annualized_last_dose_time = (
                    time - self.last_dose_time
                ) / params.model.time.steps_per_year
                annualized_half_life = params.prep.half_life / 365
                load = params.prep.peak_load * (
                    (0.5) ** (annualized_last_dose_time / annualized_half_life)
                )
                return -1.0 * np.exp(-5.528636721 * load)

        return 1.0

    # =============== HELPER METHODS ===================

    def initiate(self, model: "model.HIVModel", force: bool = False):
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

        if force:
            self.enroll(model.run_random, model.time)
        elif params.prep.target_as_prob:
            if "Racial" in params.prep.target_model:
                if (
                    model.run_random.random()
                    <= params.demographics[self.agent.race][
                        self.agent.sex_type
                    ].prep.target
                ):
                    self.enroll(model.run_random, model.time)
            else:
                if model.run_random.random() <= params.prep.target:
                    self.enroll(model.run_random, model.time)
        else:
            if "Racial" in params.prep.target_model:
                num_prep_agents = self.counts[self.agent.race]
                all_hiv_agents = model.pop.hiv_agents.members
                all_race = {
                    a for a in model.pop.all_agents if a.race == self.agent.race
                }

                num_hiv_agents = len(all_hiv_agents & all_race)
                target_prep = (len(all_race) - num_hiv_agents) * params.demographics[
                    self.agent.race
                ][self.agent.sex_type].prep.target
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
                self.enroll(model.run_random, model.time)

    def enroll(self, rand_gen, time):
        """
        Enroll an agent in PrEP

        args:
            rand_gen: random number generator
        """
        params = self.agent.location.params

        self.active = True
        self.time = time
        self.last_dose_time = time

        self.adherent = (
            rand_gen.random()
            < params.demographics[self.agent.race][self.agent.sex_type].prep.adherence
        )

        if "Inj" in params.prep.type and "Oral" in params.prep.type:
            if rand_gen.random() < params.prep.lai.prob:
                self.type = "Inj"
            else:
                self.type = "Oral"

        else:
            self.type = params.prep.type[0]

        self.add_agent(self.agent)

    def progress(self, model: "model.HIVModel", force: bool = False):
        """
        Update agent's PrEP status and discontinue stochastically or if `force` is True

        args:
            model: instance of the HIVModel being run
            force: whether to force discontinuation of PrEP
        """
        if force:
            self.discontinue()  # TO_REVIEW should this just remove the agent from counts, or discontinue? does it depend on type?
            return

        if self.type == "Oral":
            if (
                model.run_random.random()
                < self.agent.location.params.demographics[self.agent.race][
                    self.agent.sex_type
                ].prep.discontinue
            ):
                self.discontinue()
            else:
                self.last_dose_time = model.time

        # TO_REVIEW should inj prep have a way to continue at the year mark (besides maybe getting prep again through the normal channels of enrollment)?
        if (
            self.type == "Inj"
            and self.last_dose_time
            + self.agent.location.params.model.time.steps_per_year
            == model.time
        ):
            self.discontinue()

    def discontinue(self):
        """
        Discontinue PrEP usage
        """
        self.active = False
        self.type = ""
        self.time = None
        self.last_dose_time = None

        self.remove_agent(self.agent)

    def eligible(self) -> bool:
        """
        Determine if an agent is eligible for PrEP

        returns:
            whether the agent is eligible
        """
        target_model = self.agent.location.params.prep.target_model
        gender = self.agent.location.params.classes.sex_types[
            self.agent.sex_type
        ].gender

        if (
            self.active
            or self.agent.vaccine.active  # type: ignore[attr-defined]
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
            if self.agent.syringe_services.active and self.cdc_eligible():  # type: ignore[attr-defined]
                return True

        if "ssp" in target_model:
            if self.agent.syringe_services.active:  # type: ignore[attr-defined]
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
