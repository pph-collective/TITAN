# mypy: always-true=HighRisk

from typing import Dict, Set, ClassVar

from .base_feature import BaseFeature
from .. import utils
from ..agent import Agent
from ..population import Population
from ..model import HIVModel


class HighRisk(BaseFeature):

    name = "high_risk"
    stats = [
        "high_risk_new",
        "high_risk_new_hiv",
        "high_risk_new_aids",
        "high_risk_new_dx",
        "high_risk_new_haart",
        "inf_HR6m",
        "inf_HRever",
    ]
    """
        High Risk collects the following stats:

        * high_risk_new - number of agents that became active high risk this time step
        * high_risk_new_hiv - number of agents that became active high risk this time step with HIV
        * high_risk_new_aids - number of agents that became active high risk this time step with AIDS
        * high_risk_new_dx - number of agents that became active high risk this time step with diagnosed HIV
        * high_risk_new_haart - number of agents that became active high risk this time step with active HAART
        * inf_HR6m - number of agents that became active with HIV this time step who are high risk
        * inf_HRever - number of agents that became active with HIV this time step were ever high risk
    """

    new_agents: ClassVar[Set["Agent"]] = set()
    count: ClassVar[int] = 0

    def __init__(self, agent: "Agent"):
        super().__init__(agent)

        self.active = False
        self.time = 0
        self.ever = False

    def init_agent(self, pop: "Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        Based on agent demographic params, randomly initialize agent as high risk.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        if (
            pop.pop_random.random()
            < self.agent.location.params.demographics[self.agent.race][
                self.agent.sex_type
            ].high_risk.init
        ):
            self.become_high_risk()

    def update_agent(self, model: "HIVModel"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `HIVModel.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        Update high risk agents or remove them from high risk pool.  An agent becomes high_risk through the incarceration feature

        args:
            model: the instance of HIVModel currently being run
        """
        if not self.active:
            return None

        if self.time > 0:
            self.time -= 1
        else:
            self.remove_agent(self.agent)
            self.active = False

            if model.params.features.incar:
                for bond in self.agent.location.params.high_risk.partnership_types:
                    self.agent.mean_num_partners[
                        bond
                    ] -= self.agent.location.params.high_risk.partner_scale
                    self.agent.mean_num_partners[bond] = max(
                        0, self.agent.mean_num_partners[bond]
                    )  # make sure not negative
                    self.agent.target_partners[bond] = utils.poisson(
                        self.agent.mean_num_partners[bond], model.np_random
                    )
                    while (
                        len(self.agent.partners[bond])
                        > self.agent.target_partners[bond]
                    ):
                        rel = utils.safe_random_choice(
                            self.agent.relationships, model.run_random
                        )
                        if rel is not None:
                            rel.progress(force=True)
                            model.pop.remove_relationship(rel)

    @classmethod
    def add_agent(cls, agent: "Agent", new_agent: bool = True):
        """
        Add an agent to the class (not instance).

        Increment the count of high risk agents. Add the agent to the set of newly high risk agents.

        args:
            agent: the agent to add to the class attributes
            new_agent: whether the agent is newly high risk
        """
        cls.count += 1
        if new_agent:
            cls.new_agents.add(agent)

    @classmethod
    def remove_agent(cls, agent: "Agent"):
        """
        Remove an agent from the class (not instance).

        Decrement the count of high risk agents.

        args:
            agent: the agent to remove from the class attributes
        """
        cls.count -= 1

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
            stats["high_risk_new"] += 1
            if self.agent.hiv:
                stats["high_risk_new_hiv"] += 1
                if self.agent.aids:
                    stats["high_risk_new_aids"] += 1
                if self.agent.hiv_dx:
                    stats["high_risk_new_dx"] += 1
                    if self.agent.haart.active:  # type: ignore[attr-defined]
                        stats["high_risk_new_haart"] += 1

        if self.agent.hiv_time == 1:  # newly hiv
            if self.active:
                stats["inf_HR6m"] += 1
            if self.ever:
                stats["inf_HRever"] += 1

    # ============== HELPER METHODS ================

    def become_high_risk(self, duration: int = None):
        """
        Mark an agent as high risk and assign a duration to their high risk period

        args:
            duration: duration of the high risk period, defaults to param value if not passed [params.high_risk.sex_based]
        """

        if not self.agent.location.params.features.high_risk:
            return None

        is_new_agent = not self.ever
        self.add_agent(self.agent, new_agent=is_new_agent)

        self.active = True
        self.ever = True

        if duration is not None:
            self.time = duration
        else:
            self.time = self.agent.location.params.high_risk.sex_based[
                self.agent.sex_type
            ].duration
