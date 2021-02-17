from typing import Dict, ClassVar

from . import base_feature
from .. import agent
from .. import population
from .. import model
from ..parse_params import ObjMap


class HAART(base_feature.BaseFeature):
    """
    Highly Active Antiretroviral Theray (HAART) is a treatment regimen.
    """

    name = "haart"
    stats = ["haart"]
    """
        HAART collects the following stats:

        * haart - number of agents with active haart
    """

    counts: ClassVar[Dict] = {}

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)

        self.active = False
        self.ever = False
        self.adherent = False

    @classmethod
    def init_class(cls, params: "ObjMap"):
        """
        Initialize the counts dictionary for the races and sex_types in the model.

        args:
            params: the population params
        """
        cls.counts = {
            race: {sex_type: 0 for sex_type in params.classes.sex_types}
            for race in params.classes.races
        }

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        An Agent can only be initialized with HAART if they are HIV+ and diagnosed.  They are randomly assigned to HAART with a probability based on demographics and then if assigned to haart, assigned an adherence based on those same demographic params.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        haart_params = (
            self.agent.location.params.demographics[self.agent.race]
            .sex_type[self.agent.sex_type]
            .drug_type[self.agent.drug_type]
            .haart
        )
        if (
            self.agent.hiv.dx  # type: ignore[attr-defined]
            and pop.pop_random.random() < haart_params.init
        ):
            self.initiate(pop.pop_random, haart_params, "init")

    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `TITAN.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        Account for HIV treatment through highly active antiretroviral therapy (HAART). HAART was implemented in 1996, hence, there is treatment only after 1996. HIV treatment assumes that the agent knows their HIV+ status (`dx` is True).

        args:
            model: the instance of TITAN currently being run
        """
        if (
            self.agent.hiv.dx  # type: ignore[attr-defined]
            and model.time >= model.params.hiv.start_time  # haart starts with hiv
        ):
            # Determine probability of HIV treatment
            haart_params = (
                self.agent.location.params.demographics[self.agent.race]
                .sex_type[self.agent.sex_type]
                .drug_type[self.agent.drug_type]
                .haart
            )
            # Go on HAART
            if not self.active:
                self.enroll(model, haart_params)

            # Update agents on HAART
            else:
                # Go off HAART
                if model.run_random.random() < haart_params.discontinue:
                    self.active = False
                    self.adherent = False
                    self.remove_agent(self.agent)
                # Become non-adherent
                elif (
                    self.adherent
                    and model.run_random.random() < haart_params.adherence.discontinue
                ):
                    self.adherent = False
                # Become adherent
                elif (
                    not self.adherent
                    and model.run_random.random() < haart_params.adherence.become
                ):
                    self.adherent = True

    @classmethod
    def add_agent(cls, agent: "agent.Agent"):
        """
        Add an agent to the class (not instance).

        Increments `counts` or haart agents by race and sex_type for the given agent.

        args:
            agent: the agent to add to the class attributes
        """
        cls.counts[agent.race][agent.sex_type] += 1

    @classmethod
    def remove_agent(cls, agent: "agent.Agent"):
        """
        Remove an agent from the class (not instance).

        Decrements `counts` or haart agents by race and sex_type for the given agent.

        args:
            agent: the agent to remove from the class attributes
        """
        cls.counts[agent.race][agent.sex_type] -= 1

    def set_stats(self, stats: Dict[str, int], time: int):
        if self.active:
            stats["haart"] += 1

    def get_transmission_risk_multiplier(self, time: int, interaction_type: str):
        """
        Get a multiplier for how haart reduces hiv transmission risk based on interaction type and params.

        By default, returns 1.0

        args:
            time: the current model time step
            interaction_type: The type of interaction where the agent could transmit HIV (e.g. 'sex', 'injection' - from [params.classes.interaction_types])
        """
        prob = 1.0
        if self.active:
            params = self.agent.location.params
            adherence = "adherent" if self.adherent else "non_adherent"
            if interaction_type == "injection":
                prob = params.partnership.injection.transmission.haart_scaling[
                    adherence
                ]
            elif interaction_type == "sex":
                prob = params.partnership.sex.haart_scaling[self.agent.sex_type][
                    adherence
                ]

            # Tuning parameter for ART efficiency
            prob *= params.calibration.haart.transmission

        return prob

    def aids_scale(self):
        prob = 1.0
        if self.active:
            adherence = "adherent" if self.adherent else "non_adherent"
            prob = self.agent.location.params.haart.aids_scale[adherence]

        return prob

    # =========== HELPER METHODS ============

    def enroll(self, model: "model.TITAN", haart_params: ObjMap):
        """
        Determine whether to enroll an agent in HAART.

        args:
            model: the instance of TITAN currently being run
            haart_params: the HAART demographic params for this agent
        """
        if self.agent.location.params.haart.use_cap:
            self.enroll_cap(model, haart_params)
        else:
            self.enroll_prob(model, haart_params)

    def enroll_cap(self, model: "model.TITAN", haart_params: ObjMap):
        """
        Determine whether to enroll an agent in HAART using the cap method.

        args:
            model: the instance of TITAN currently being run
            haart_params: the HAART demographic params for this agent
        """
        race = self.agent.race
        sex_type = self.agent.sex_type
        # if HAART is based on cap instead of prob, determine number of
        # HAART agents based on % of diagnosed agents
        num_dx_agents = self.agent.hiv.dx_counts[race][sex_type]  # type: ignore[attr-defined]
        num_haart_agents = self.counts[race][sex_type]

        # take value from dictionary for cap
        if num_haart_agents < (haart_params.cap * num_dx_agents):
            self.initiate(model.run_random, haart_params, "prob")

    def enroll_prob(self, model: "model.TITAN", haart_params: ObjMap):
        """
        Determine whether to enroll an agent in HAART using probability method.

        args:
            model: the instance of TITAN currently being run
            haart_params: the HAART demographic params for this agent
        """
        if self.ever and self.agent.location.params.haart.use_reinit:
            if model.run_random.random() < haart_params.reinit.prob:
                self.initiate(model.run_random, haart_params, "prob")
        else:
            # Find enroll probability based on time since diagnosis
            enroll_prob = 0.0
            dx_duration = model.time - self.agent.hiv.dx_time  # type: ignore[attr-defined]
            for i in haart_params.enroll.values():
                if i.start <= dx_duration < i.stop:
                    enroll_prob = i.prob * model.calibration.haart.coverage
                    break

            if model.run_random.random() < (enroll_prob):
                self.initiate(model.run_random, haart_params, "prob")

    def initiate(self, rand_gen, haart_params: ObjMap, init_or_prob: str):
        """
        Initiate an agent with HAART and add them to the population.

        args:
            model: the instance of TITAN currently being run
        """
        self.adherent = rand_gen.random() < haart_params.adherence[init_or_prob]

        # Add agent to HAART class set, update agent params
        self.active = True
        self.ever = True
        self.add_agent(self.agent)
