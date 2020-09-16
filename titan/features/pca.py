from .base_feature import BaseFeature


class PCA(BaseFeature):

    name = "pca"
    stats = []

    def __init__(self, agent):
        super().__init__(agent)

        self.active = False
        self.suitable = False
        self.awareness = False
        self.opinion = 0.0

    def update_agent(self, model):
        pca_params = self.agent.location.params.pca
        if (
            model.time >= pca_params.start_time
            and model.run_random.random() < pca_params.awareness.prob
        ):
            self.awareness = True
            if model.run_random.random() < pca_params.prep.prob:
                self.agent.prep.initiate(model, force=True)

    def init_agent(self, pop, time):
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
