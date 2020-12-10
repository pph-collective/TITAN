from typing import Dict, Optional

from . import base_feature
from .. import agent
from .. import population
from .. import model


class Incar(base_feature.BaseFeature):

    name = "incar"
    stats = ["incar", "incar_hiv", "new_release", "new_release_hiv"]
    """
        Incar collects the following stats:

        * incar - number of agents with active incar
        * incar_hiv - number of agents with active incar and HIV
        * new_release - number of agents released this timestep
        * new_release_hiv - number of agents releasted this timestep with HIV
    """

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)

        self.active = False
        self.time: Optional[int] = None
        self.release_time: Optional[int] = None

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        Run incarceration assignment on an agent.  The duration of incarceration at initialization is different than the ongoing to reflect that agents with longer durations will be more highly represented in that population at any given point in time.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        incar_params = (
            self.agent.location.params.demographics[self.agent.race]
            .sex_type[self.agent.sex_type]
            .incar
        )
        jail_duration = incar_params.duration.init

        prob_incar = incar_params.init
        if pop.pop_random.random() < prob_incar:
            self.active = True
            bin = 1
            current_p_value = jail_duration[bin].prob
            p = pop.pop_random.random()

            while p > current_p_value:
                bin += 1
                current_p_value += jail_duration[bin].prob

            self.time = time
            self.release_time = time + pop.pop_random.randrange(
                jail_duration[bin].min, jail_duration[bin].max
            )

    def update_agent(self, model: "model.HIVModel"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `HIVModel.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        Incarcerate an agent or update their incarceration variables

        args:
            model: the instance of HIVModel currently being run
        """
        hiv_bool = self.agent.hiv

        if hiv_bool:
            hiv_multiplier = self.agent.location.params.incar.hiv.multiplier
        else:
            hiv_multiplier = 1.0

        # agent is incarcerated
        if self.active:

            # Release agent
            if self.release_time == model.time:
                self.active = False

                # does agent stay on haart
                if hiv_bool:
                    if self.agent.haart.active:  # type: ignore[attr-defined]
                        if (
                            model.run_random.random()
                            <= self.agent.location.params.incar.haart.discontinue
                        ):
                            self.agent.haart.active = False  # type: ignore[attr-defined]
                            self.agent.haart.adherent = False  # type: ignore[attr-defined]

        # should the agent become incarcerated?
        elif model.run_random.random() < (
            self.agent.location.params.demographics[self.agent.race]
            .sex_type[self.agent.sex_type]
            .incar.prob
            * hiv_multiplier
            * model.calibration.incarceration
        ):
            incar_duration = (
                self.agent.location.params.demographics[self.agent.race]
                .sex_type[self.agent.sex_type]
                .incar.duration.prob
            )

            bin = current_p_value = 1
            p = model.run_random.random()
            while p >= current_p_value:
                current_p_value += incar_duration[bin].prob
                bin += 1

            self.time = model.time
            self.release_time = model.time + model.run_random.randint(
                incar_duration[bin].min, incar_duration[bin].max
            )
            self.active = True

            if hiv_bool:
                if not self.agent.hiv_dx:
                    if (
                        model.run_random.random()
                        < self.agent.location.params.incar.hiv.dx
                    ):
                        self.agent.hiv_dx = True
                else:  # Then tested and HIV, check to enroll in ART
                    if (
                        model.run_random.random()
                        < self.agent.location.params.incar.haart.prob
                    ):
                        self.agent.haart.adherent = model.run_random.random() < self.agent.location.params.incar.haart.adherence  # type: ignore[attr-defined]
                        # Add agent to HAART class set, update agent params
                        self.agent.haart.active = True  # type: ignore[attr-defined]

    def set_stats(self, stats: Dict[str, int], time: int):
        if self.release_time == time:
            stats["new_release"] += 1
            if self.agent.hiv:
                stats["new_release_hiv"] += 1

        if self.active:
            stats["incar"] += 1
            if self.agent.hiv:
                stats["incar_hiv"] += 1

    # ============== HELPER METHODS ================
