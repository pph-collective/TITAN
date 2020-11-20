from . import base_feature
from .. import agent
from .. import population
from .. import model


class ExternalExposure(base_feature.BaseFeature):

    name = "external_exposure"

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)

        self.active = False

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        If an agent has defined sex_type, with a random probability from params, assign them to be an agent with external exposure.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        params = self.agent.location.params.external_exposure
        if self.agent.sex_type == params.sex_type:
            if pop.pop_random.random() < params.init:
                self.active = True

    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `TITAN.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        If the agent has external exposure, with a probability from params, convert the agent.

        args:
            model: the instance of TITAN currently being run
        """
        params = self.agent.location.params.external_exposure
        if self.active and model.run_random.random() < params.convert_prob:
            agent_exposure = getattr(self.agent, params.exposure)
            agent_exposure.convert(model)
