from .base_feature import BaseFeature


class HAART(BaseFeature):

    name = "haart"
    stats = ["numART"]

    counts = None

    def __init__(self, agent):
        super().__init__(agent)

        self.active = False
        self.time = 0
        self.adherence = 0

    def init_agent(self, pop, time):
        agent_params = self.agent.location.params.demographics[self.agent.race][
            self.agent.population
        ]
        if (
            self.agent.hiv
            and self.agent.hiv_dx
            and pop.pop_random.random() < agent_params.haart.init
        ):
            self.active = True
            self.time = 0
            self.add_agent(self.agent)

            haart_adh = agent_params.haart.adherence
            if pop.pop_random.random() < haart_adh:
                self.adherence = 5
            else:
                self.adherence = pop.pop_random.randint(1, 4)

    def update_agent(self, model):
        if self.agent.hiv and model.time >= model.params.hiv.start_time:
            self.update_haart(model)

    @classmethod
    def add_agent(cls, agent):
        # set up if this is the first time being called
        if cls.counts is None:
            cls.init_class(agent.location.params)

        cls.counts[agent.race][agent.sex_type] += 1

    def remove_agent(cls, agent):
        cls.counts[agent.race][agent.sex_type] -= 1

    def set_stats(self, stats):
        if self.active:
            stats["numART"] += 1

    @classmethod
    def init_class(cls, params):
        cls.counts = {
            race: {sex_type: 0 for sex_type in params.classes.sex_types}
            for race in params.classes.races
        }

    def initiate(self, model):
        haart_adh = self.agent.location.params.demographics[self.agent.race][
            self.agent.sex_type
        ].haart.adherence
        if model.run_random.random() < haart_adh:
            adherence = 5
        else:
            adherence = model.run_random.randint(1, 4)

        # Add agent to HAART class set, update agent params
        self.active = True
        self.adherence = adherence
        self.time = model.time
        self.add_agent(self.agent)

    def update_haart(self, model):
        """
        Account for HIV treatment through highly active antiretroviral therapy
            (HAART).
            HAART was implemented in 1996, hence, there is treatment only after 1996.
            HIV treatment assumes that the agent knows their HIV+ status (`dx` is True).

        args:
            agent: agent being updated
        """
        # Check valid input
        assert self.agent.hiv

        if self.counts is None:
            self.init_class(model.params)

        # Determine probability of HIV treatment
        if self.agent.hiv_dx:
            haart_params = self.agent.location.params.demographics[self.agent.race][
                self.agent.sex_type
            ].haart
            # Go on HAART
            if not self.active:
                if self.agent.location.params.hiv.haart_cap:
                    # if HAART is based on cap instead of prob, determine number of
                    # HAART agents based on % of diagnosed agents
                    num_dx_agents = model.pop.dx_counts[self.agent.race][
                        self.agent.sex_type
                    ]
                    num_haart_agents = self.counts[self.agent.race][self.agent.sex_type]

                    if num_haart_agents < (haart_params.prob * num_dx_agents):
                        self.initiate(model)
                else:
                    if model.run_random.random() < (
                        haart_params.prob * model.calibration.haart.coverage
                    ):
                        self.initiate(model)
            # Go off HAART
            elif self.active and model.run_random.random() < haart_params.discontinue:
                self.active = False
                self.adherence = 0
                self.time = 0
                self.remove_agent(self.agent)
