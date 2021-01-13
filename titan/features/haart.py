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
        agent_params = (
            self.agent.location.params.demographics[self.agent.race]
            .sex_type[self.agent.sex_type]
            .drug_type[self.agent.drug_type]
        )
        if (
            self.agent.hiv.active  # type: ignore[attr-defined]
            and self.agent.hiv.dx  # type: ignore[attr-defined]
            and pop.pop_random.random() < agent_params.haart.init
        ):
            self.active = True
            self.ever = True
            self.add_agent(self.agent)

            haart_adh = agent_params.haart.adherence.init
            if pop.pop_random.random() < haart_adh:
                self.adherent = True

    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `TITAN.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        Account for HIV treatment through highly active antiretroviral therapy (HAART). HAART was implemented in 1996, hence, there is treatment only after 1996. HIV treatment assumes that the agent knows their HIV+ status (`dx` is True).

        args:
            model: the instance of TITAN currently being run
        """
        if (
            self.agent.hiv.active  # type: ignore[attr-defined]
            and self.agent.hiv.dx  # type: ignore[attr-defined]
            and model.time >= model.params.hiv.start_time
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
                if self.agent.location.params.hiv.haart_cap:
                    # if HAART is based on cap instead of prob, determine number of
                    # HAART agents based on % of diagnosed agents
                    num_dx_agents = self.agent.hiv.dx_counts[self.agent.race][  # type: ignore[attr-defined]
                        self.agent.sex_type
                    ]
                    num_haart_agents = self.counts[self.agent.race][self.agent.sex_type]

                    # take value from dictionary for cap
                    if num_haart_agents < (haart_params.cap * num_dx_agents):
                        self.initiate(model)
                else:
                    if self.ever and self.agent.location.params.hiv.reinit_haart:
                        if model.run_random.random() < haart_params.reinit.prob:
                            self.initiate(model)
                    else:
                        enroll_prob = 0
                        # Find enroll probability based on time since diagnosis
                        for i in haart_params.prob.values():
                            if i.start <= (model.time - self.agent.hiv.dx_time) < i.stop:  # type: ignore[attr-defined]
                                enroll_prob = i.enroll_prob
                                break

                        if model.run_random.random() < (
                            enroll_prob * model.calibration.haart.coverage
                        ):
                            self.initiate(model)

            # Go off HAART
            elif self.active:
                if model.run_random.random() < haart_params.discontinue:
                    self.active = False
                    self.adherent = False
                    self.remove_agent(self.agent)
                elif (
                    self.adherent
                    and model.run_random.random() < haart_params.adherence.discontinue
                ):
                    self.adherent = False
                elif (
                    not self.adherent
                    and model.run_random.random() < haart_params.adherence.prob
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

    # =========== HELPER METHODS ============

    def initiate(self, model: "model.TITAN"):
        """
        Initiate an agent with HAART and add them to the population.

        args:
            model: the instance of TITAN currently being run
        """
        self.adherent = (
            model.run_random.random()
            < self.agent.location.params.demographics[self.agent.race]
            .sex_type[self.agent.sex_type]
            .drug_type[self.agent.drug_type]
            .haart.adherence.prob
        )

        # Add agent to HAART class set, update agent params
        self.active = True
        self.ever = True
        self.add_agent(self.agent)
