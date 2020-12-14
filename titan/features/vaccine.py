from typing import Dict, Optional

import numpy as np  # type: ignore

from . import base_feature
from .. import agent
from .. import population
from .. import model


class Vaccine(base_feature.BaseFeature):

    name = "vaccine"
    stats = ["vaccine"]
    """
        Vaccine collects the following stats:

        * vaccine - number of agents with active vaccine
    """

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)
        self.active = False
        self.time: Optional[int] = None
        self.type = ""

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        If the agent is HIV-, randomly vaccinate per the params.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        if (
            not self.agent.hiv.active  # type: ignore[attr-defined]
            and self.agent.location.params.vaccine.on_init
            and pop.pop_random.random()
            < self.agent.location.params.demographics[self.agent.race]
            .sex_type[self.agent.sex_type]
            .vaccine.init
        ):
            self.vaccinate(time)

    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `TITAN.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        If PrEP feature is enable and the agent is not active PrEP and not HIV, either update or stochastically vaccinate the agent.

        args:
            model: the instance of TITAN currently being run
        """
        if (
            not self.agent.prep.active  # type: ignore[attr-defined]
            and not self.agent.hiv.active  # type: ignore[attr-defined]
        ):
            vaccine_params = self.agent.location.params.vaccine
            agent_params = (
                self.agent.location.params.demographics[self.agent.race]
                .sex_type[self.agent.sex_type]
                .vaccine
            )

            if self.active:
                if (
                    vaccine_params.booster
                    and (model.time - self.time) == agent_params.booster.interval
                    and model.run_random.random() < agent_params.booster.prob
                ):
                    self.vaccinate(model.time)
            elif model.time == vaccine_params.start_time:
                if model.run_random.random() < agent_params.prob:
                    self.vaccinate(model.time)

    def set_stats(self, stats: Dict[str, int], time: int):
        if self.active:
            stats["vaccine"] += 1

    def get_acquisition_risk_multiplier(self, time: int, interaction_type: str):
        """
        Get a multiplier for how vaccine affects acquisition of HIV for the given interaction_type.

        By default, returns 1.0

        args:
            time: the current model time step
            interaction_type: The type of interaction where the agent could acquire HIV (e.g. 'sex', 'injection' - from [params.classes.interaction_types])
        """
        # not protected the time step the agent is vaccinaetd
        if self.active and self.time is not None and self.time < time:
            vaccine_time_months = (
                (time - self.time)
                / self.agent.location.params.model.time.steps_per_year
            ) * 12

            if self.type == "HVTN702":
                return np.exp(
                    -2.88 + 0.76 * (np.log((vaccine_time_months + 0.001) * 30))
                )
            elif self.type == "RV144":
                return np.exp(-2.40 + 0.76 * (np.log(vaccine_time_months)))

        return 1.0

    # ============= HELPER METHODS =============

    def vaccinate(self, time):
        """
        Vaccinate an agent and update relevant fields.
        """
        self.active = True
        self.type = self.agent.location.params.vaccine.type
        self.time = time
