class BaseFeature:
    def __init__(self):
        self.active = False

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
