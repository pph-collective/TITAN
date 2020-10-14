from typing import Dict, ClassVar, Set

from . import base_feature
from .. import utils
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

    new_releases: ClassVar[Set["agent.Agent"]] = set()

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)

        self.active = False
        self.duration = 0

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        Run incarceration assignment on an agent.  The duration of incarceration at initialization is different than the ongoing to reflect that agents with longer durations will be more highly represented in that population at any given point in time.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        incar_params = self.agent.location.params.demographics[self.agent.race][
            self.agent.sex_type
        ].incar
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

            self.duration = pop.pop_random.randrange(
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
            self.duration -= 1

            # Release agent
            if self.duration == 0:
                self.add_agent(self.agent)
                self.active = False

                # become high risk on release
                if (
                    not self.agent.high_risk.active and model.params.features.high_risk  # type: ignore[attr-defined]
                ):  # If behavioral treatment on and agent HIV, ignore HR period.
                    self.agent.high_risk.become_high_risk(model.time)  # type: ignore[attr-defined]
                    for bond in self.agent.location.params.high_risk.partnership_types:
                        self.agent.mean_num_partners[
                            bond
                        ] += self.agent.location.params.high_risk.partner_scale
                        self.agent.target_partners[bond] = utils.poisson(
                            self.agent.mean_num_partners[bond], model.np_random
                        )
                    model.pop.update_partnerability(self.agent)

                # does agent stay on haart
                if hiv_bool:
                    if self.agent.haart.active:  # type: ignore[attr-defined]
                        if (
                            model.run_random.random()
                            <= self.agent.location.params.incar.haart.discontinue
                        ):
                            self.agent.haart.active = False  # type: ignore[attr-defined]
                            self.agent.haart.adherence = 0  # type: ignore[attr-defined]

                        # END FORCE

        # should the agent become incarcerated?
        elif model.run_random.random() < (
            self.agent.location.params.demographics[self.agent.race][
                self.agent.sex_type
            ].incar.prob
            * hiv_multiplier
            * model.calibration.incarceration
        ):
            incar_duration = self.agent.location.params.demographics[self.agent.race][
                self.agent.sex_type
            ].incar.duration.prob

            bin = current_p_value = 1
            p = model.run_random.random()
            while p >= current_p_value:
                current_p_value += incar_duration[bin].prob
                bin += 1

            self.duration = model.run_random.randint(
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
                        if (
                            model.run_random.random()
                            < self.agent.location.params.incar.haart.adherence
                        ):
                            adherence = 5
                        else:
                            adherence = model.run_random.randint(1, 4)

                        # Add agent to HAART class set, update agent params
                        self.agent.haart.active = True  # type: ignore[attr-defined]
                        self.agent.haart.adherence = adherence  # type: ignore[attr-defined]

            # PUT PARTNERS IN HIGH RISK
            if model.params.features.high_risk:
                for partner in self.agent.get_partners(
                    self.agent.location.params.high_risk.partnership_types
                ):
                    if not partner.high_risk.active:  # type: ignore[attr-defined]
                        if (
                            model.run_random.random()
                            < partner.location.params.high_risk.prob
                        ):
                            partner.high_risk.become_high_risk(model.time)  # type: ignore[attr-defined]

    @classmethod
    def add_agent(cls, agent: "agent.Agent"):
        cls.new_releases.add(agent)

    @classmethod
    def update_pop(cls, model: "model.HIVModel"):
        """
        Update the feature for the entire population (class method).

        This is called in `HIVModel.update_all_agents` before agent-level updates are made.

        Resets the tracking set for `new_releases` as this is called before `update_agent`

        args:
            model: the instance of HIVModel currently being run
        """
        # population is updated before agents, so clear set at the beginning of updates
        cls.new_releases = set()

    def set_stats(self, stats: Dict[str, int], time: int):
        if self.agent in self.new_releases:
            stats["new_release"] += 1
            if self.agent.hiv:
                stats["new_release_hiv"] += 1

        if self.active:
            stats["incar"] += 1
            if self.agent.hiv:
                stats["incar_hiv"] += 1

    # ============== HELPER METHODS ================
