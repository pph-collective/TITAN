from .base_feature import BaseFeature


class MSMW(BaseFeature):

    name = "msmw"

    def __init__(self, agent: "Agent"):
        super().__init__(agent)

        self.active = False

    def init_agent(self, pop: "Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        If an agent has `sex_type == "HM"`, with a random probability from params, assign them to be a Man who Sleeps with Men and Women (MSMW).

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        if self.agent.sex_type == "HM":
            if pop.pop_random.random() < self.agent.location.params.msmw.prob:
                self.active = True

    def update_agent(self, model: "HIVModel"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `HIVModel.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        If the agent is MSMW, with a probability from params, hiv convert the agent.

        args:
            model: the instance of HIVModel currently being run
        """
        if (
            self.active
            and model.run_random.random() < self.agent.location.params.msmw.hiv.prob
        ):
            model.hiv_convert(self.agent)
