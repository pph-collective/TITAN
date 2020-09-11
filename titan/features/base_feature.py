class BaseFeature:

    stats = []
    name = ""

    def __init__(self, agent):
        self.active = False
        self.agent = agent

    def update_agent(self, model):
        pass

    def init_agent(self, pop, time):
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

    def set_stats(self, stats):
        pass
