from .base_feature import BaseFeature

class HighRisk(BaseFeature):

    new_agents = set()
    count = 0

    def __init__(self, agent):
        super().__init__(agent)

        self.active = False
        self.time = 0
        self.ever = False

    @classmethod
    def update_pop(cls, model):
        # population is updated before agents, so clear set at the beginning of updates
        cls.new_agents = set()

    @classmethod
    def add_agent(cls, agent, new_agent: bool = True):
        cls.count += 1
        if new_agent:
            cls.new_agents.add(agent)

    @classmethod
    def remove_agent(cls, agent):
        cls.count -= 1

    def update_agent(self, model):
        """
        Update high risk agents or remove them from high risk pool.  An agent
        becomes high_risk through the incarceration feature
        """
        if not self.active:
            return None

        if self.time > 0:
            self.time -= 1
        else:
            self.remove_agent(agent)
            self.high_risk = False

            if model.features.incar:
                for bond in self.agent.location.params.high_risk.partnership_types:
                    agent.mean_num_partners[
                        bond
                    ] -= self.agent.location.params.high_risk.partner_scale
                    self.agent.mean_num_partners[bond] = max(
                        0, self.agent.mean_num_partners[bond]
                    )  # make sure not negative
                    self.agent.target_partners[bond] = utils.poisson(
                        self.agent.mean_num_partners[bond], model.np_random
                    )
                    while len(self.agent.partners[bond]) > self.agent.target_partners[bond]:
                        rel = utils.safe_random_choice(
                            self.agent.relationships, model.run_random
                        )
                        if rel is not None:
                            rel.progress(force=True)
                            model.pop.remove_relationship(rel)

    def become_high_risk(self, model, duration: int = None):
        """
        Mark an agent as high risk and assign a duration to their high risk period

        args:
            agent: agent becoming high risk
            duration: duration of the high risk period, defaults to param value if not passed [params.high_risk.sex_based]
        """

        if not self.features.high_risk:
            return None

        is_new_agent = not self.ever
        self.add_agent(agent, new_agent=is_new_agent)

        self.active = True
        self.ever = True

        if duration is not None:
            self.time = duration
        else:
            self.time = self.agent.location.params.high_risk.sex_based[
                agent.sex_type
            ].duration
