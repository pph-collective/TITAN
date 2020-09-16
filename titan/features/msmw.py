from .base_feature import BaseFeature


class MSMW(BaseFeature):

    name = "msmw"

    def __init__(self, agent):
        super().__init__(agent)

        self.active = False

    def init_agent(self, pop, time):
        if self.agent.sex_type == "HM":
            if pop.pop_random.random() < self.agent.location.params.msmw.prob:
                self.active = True

    def update_agent(self, model):
        if (
            self.active
            and model.run_random.random() < self.agent.location.params.msmw.hiv.prob
        ):
            model.hiv_convert(self.agent)
