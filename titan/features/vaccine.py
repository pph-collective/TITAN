from .base_feature import BaseFeature


class Vaccine(BaseFeature):
    def __init__(self):
        self.active = False
        self.time = 0
        self.type = ""

    def update_agent(self, agent, model):
        if model.features.prep and not agent.prep.active:
            self.advance_vaccine(agent, model)

    def advance_vaccine(self, agent, model):
        """
        Progress vaccine. Agents may receive injection or progress in time
            since injection.

        args:
            agent: agent being updated
        """
        vaccine_params = agent.location.params.vaccine
        agent_params = agent.location.params.demographics[agent.race][
            agent.sex_type
        ].vaccine

        if self.active:
            self.time += 1
            if (
                vaccine_params.booster
                and self.time == agent_params.booster.interval
                and model.run_random.random() < agent_params.booster.prob
            ):
                self.vaccinate(agent)
        elif model.time == vaccine_params.start_time:
            if model.run_random.random() < agent_params.prob:
                self.vaccinate(agent)

    def vaccinate(self, agent):
        """
        Vaccinate an agent and update relevant fields.
        """
        self.active = True
        self.type = agent.location.params.vaccine.type
        self.time = 1
