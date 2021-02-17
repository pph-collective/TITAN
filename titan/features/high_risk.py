# mypy: always-true=HighRisk

from typing import Dict, Optional

from . import base_feature
from .. import utils
from .. import agent
from .. import population
from .. import model


class HighRisk(base_feature.BaseFeature):

    name = "high_risk"
    stats = [
        "high_risk_new",
        "high_risk_new_hiv",
        "high_risk_new_aids",
        "high_risk_new_dx",
        "high_risk_new_haart",
        "hiv_new_high_risk",
        "hiv_new_high_risk_ever",
    ]
    """
        High Risk collects the following stats:

        * high_risk_new - number of agents that became active high risk this time step
        * high_risk_new_hiv - number of agents that became active high risk this time step with HIV
        * high_risk_new_aids - number of agents that became active high risk this time step with AIDS
        * high_risk_new_dx - number of agents that became active high risk this time step with diagnosed HIV
        * high_risk_new_haart - number of agents that became active high risk this time step with active HAART
        * hiv_new_high_risk - number of agents that became active with HIV this time step who are high risk
        * hiv_new_high_risk_ever - number of agents that became active with HIV this time step were ever high risk
    """

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)

        self.active = False
        self.time: Optional[int] = None
        self.duration = 0
        self.ever = False

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        Based on agent demographic params, randomly initialize agent as high risk.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        if (
            pop.pop_random.random()
            < self.agent.location.params.demographics[self.agent.race]
            .sex_type[self.agent.sex_type]
            .high_risk.init
        ):
            self.become_high_risk(pop, time)

    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `TITAN.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        Update high risk agents or remove them from high risk pool.  An agent becomes high_risk through the incarceration feature

        args:
            model: the instance of TITAN currently being run
        """
        if not self.active:
            # released last step, evaluate agent for high risk
            if self.agent.incar.release_time == model.time - 1:  # type: ignore[attr-defined]
                self.become_high_risk(model.pop, model.time)

            # incarcerated last step, evaluate agent's partners for high risk
            elif self.agent.incar.time == model.time - 1:  # type: ignore[attr-defined]

                # put partners in high risk
                for partner in self.agent.get_partners(
                    self.agent.location.params.high_risk.partnership_types
                ):
                    if (
                        not partner.high_risk.active  # type: ignore[attr-defined]
                        and model.run_random.random()
                        < partner.location.params.high_risk.prob
                    ):
                        partner.high_risk.become_high_risk(model.pop, model.time)  # type: ignore[attr-defined]
        elif self.duration > 0:
            self.duration -= 1
        else:
            self.active = False

            self.update_partner_numbers(
                model.pop, -1 * self.agent.location.params.high_risk.partner_scale
            )

            for bond in self.agent.location.params.high_risk.partnership_types:
                num_ended = 0
                while (
                    len(self.agent.partners[bond]) - num_ended
                ) > self.agent.target_partners[bond]:
                    rel = utils.safe_random_choice(
                        self.agent.relationships, model.run_random
                    )
                    if rel is not None:
                        num_ended += 1
                        rel.duration = 0  # will end on next step

    def set_stats(self, stats: Dict[str, int], time: int):
        if self.time == time:
            stats["high_risk_new"] += 1
            if self.agent.hiv.active:  # type: ignore[attr-defined]
                stats["high_risk_new_hiv"] += 1
                if self.agent.hiv.aids:  # type: ignore[attr-defined]
                    stats["high_risk_new_aids"] += 1
                if self.agent.hiv.dx:  # type: ignore[attr-defined]
                    stats["high_risk_new_dx"] += 1
                    if self.agent.haart.active:  # type: ignore[attr-defined]
                        stats["high_risk_new_haart"] += 1

        # newly hiv
        if self.agent.hiv.time == time:  # type: ignore[attr-defined]
            if self.active:
                stats["hiv_new_high_risk"] += 1
            if self.ever:
                stats["hiv_new_high_risk_ever"] += 1

    # ============== HELPER METHODS ================

    def become_high_risk(
        self, pop: "population.Population", time: int, duration: int = None
    ):
        """
        Mark an agent as high risk and assign a duration to their high risk period

        args:
            pop: the model poopulation
            time: the time step the agent is becoming high risk
            duration: duration of the high risk period, defaults to param value if not passed [params.high_risk.sex_based]
        """

        if not self.agent.location.params.features.high_risk:
            return None

        if not self.ever:
            self.time = time

        self.active = True
        self.ever = True

        if duration is not None:
            self.duration = duration
        else:
            self.duration = self.agent.location.params.high_risk.sex_based[
                self.agent.sex_type
            ].duration

        self.update_partner_numbers(
            pop, self.agent.location.params.high_risk.partner_scale
        )

    def update_partner_numbers(self, pop: "population.Population", amount: int):
        """
        Update the agent's mean and target partner numbers by the amount passed.  Update partnerability for the population.

        args:
            pop: the model population
            amount: the positive or negatative amount to adjust the mean by
        """
        for bond in self.agent.location.params.high_risk.partnership_types:
            self.agent.mean_num_partners[bond] += amount  # could be negative
            self.agent.target_partners[bond] = utils.poisson(
                self.agent.mean_num_partners[bond], pop.np_random
            )
            pop.update_partnerability(self.agent)
