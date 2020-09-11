from .base_feature import BaseFeature


class Vaccine(BaseFeature):

    stats = ["vaccinated"]
    name = "vaccine"

    def __init__(self, agent):
        super().__init__(agent)
        self.active = False
        self.time = 0
        self.type = ""

    def init_agent(self, pop, time):
        if (
            not self.agent.hiv
            and self.agent.location.params.vaccine.on_init
            and pop.pop_random.random()
            < self.agent.location.params.demographics[self.agent.race][
                self.agent.sex_type
            ].vaccine.init
        ):
            self.vaccinate()

    def update_agent(self, model):
        if (
            model.params.features.prep
            and not self.agent.prep.active
            and not self.agent.hiv
        ):
            self.advance_vaccine(model)

    def set_stats(self, stats):
        if self.active:
            stats["vaccinated"] += 1

    def advance_vaccine(self, model):
        """
        Progress vaccine. Agents may receive injection or progress in time
            since injection.

        args:
            agent: agent being updated
        """
        vaccine_params = self.agent.location.params.vaccine
        agent_params = self.agent.location.params.demographics[self.agent.race][
            self.agent.sex_type
        ].vaccine

        if self.active:
            self.time += 1
            if (
                vaccine_params.booster
                and self.time == agent_params.booster.interval
                and model.run_random.random() < agent_params.booster.prob
            ):
                self.vaccinate()
        elif model.time == vaccine_params.start_time:
            if model.run_random.random() < agent_params.prob:
                self.vaccinate()

    def vaccinate(self):
        """
        Vaccinate an agent and update relevant fields.
        """
        self.active = True
        self.type = self.agent.location.params.vaccine.type
        self.time = 1
