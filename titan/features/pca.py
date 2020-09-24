from .base_feature import BaseFeature
from ..agent import Agent
from ..population import Population
from ..model import HIVModel


class PCA(BaseFeature):

    name = "pca"

    def __init__(self, agent: "Agent"):
        super().__init__(agent)

        self.active = False
        self.suitable = False
        self.awareness = False
        self.opinion = 0.0

    def init_agent(self, pop: "Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        With a random probability from `params.pca`, make agent aware and assign an opinion.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        pca_params = self.agent.location.params.pca
        if pop.pop_random.random() < pca_params.awareness.init:
            self.awareness = True

            attprob = pop.pop_random.random()
            pvalue = 0.0
            for bin, fields in pca_params.attitude.items():
                pvalue += fields.prob
                if attprob < pvalue:
                    self.opinion = bin
                break

    def update_agent(self, model: "HIVModel"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `HIVModel.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        If the model's time step is at least the PCA start time, randomly make agent aware and then randomly initiate PrEP on the agent.

        args:
            model: the instance of HIVModel currently being run
        """
        pca_params = self.agent.location.params.pca
        if (
            model.time >= pca_params.start_time
            and model.run_random.random() < pca_params.awareness.prob
        ):
            self.awareness = True
            if model.run_random.random() < pca_params.prep.prob:
                self.agent.prep.initiate(model, force=True)  # type: ignore[attr-defined]
