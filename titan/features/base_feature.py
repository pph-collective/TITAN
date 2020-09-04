class BaseFeature:

    stats = []

    def __init__(self, agent):
        self.active = False
        self.agent = agent

    def update_agent(self, agent, model):
        pass

    def init_agent(self, agent, model):
        pass

    @classmethod
    def add_agent(cls, agent):
        pass

    @classmethod
    def remove_agent(cls, agent):
        pass

    @classmethod
    def update_pop(cls, model):
        pass

    def set_stats(self, stats, agent):
        pass
