from .base_feature import BaseFeature
from .. import utils


class Incar(BaseFeature):

    name = "incar"
    stats = ["incar", "incarHIV", "newRelease", "newReleaseHIV"]

    new_releases = set()

    def __init__(self, agent):
        super().__init__(agent)

        self.active = False
        self.time = 0

    def update_agent(self, model):
        self.incarcerate(model)

    def init_agent(self, pop, time):
        """
        Run incarceration assignment on a population.  The duration of incarceration at initialization is different than the ongoing to reflect that agents with longer durations will be more highly represented in that population at any given point in time.
        """
        incar_params = self.agent.location.params.demographics[self.agent.race][
            self.agent.sex_type
        ].incar
        jail_duration = incar_params.duration.init

        prob_incar = incar_params.init
        if pop.pop_random.random() < prob_incar:
            self.active = True
            bin = current_p_value = 0
            p = pop.pop_random.random()

            while p > current_p_value:
                bin += 1
                current_p_value += jail_duration[bin].prob

            self.time = pop.pop_random.randrange(
                jail_duration[bin].min, jail_duration[bin].max
            )

    @classmethod
    def update_pop(cls, model):
        # population is updated before agents, so clear set at the beginning of updates
        cls.new_releases = set()

    @classmethod
    def add_new_release_agent(cls, agent):
        cls.new_releases.add(agent)

    def set_stats(self, stats):
        if self.agent in self.new_releases:
            stats["newRelease"] += 1
            if self.agent.hiv:
                stats["newReleaseHIV"] += 1

        if self.active:
            stats["incar"] += 1
            if self.agent.hiv:
                stats["incarHIV"] += 1

    def incarcerate(self, model):
        """
        Incarcerate an agent or update their incarceration variables

        args:
            agent: agent being updated
        """
        hiv_bool = self.agent.hiv

        if hiv_bool:
            hiv_multiplier = self.agent.location.params.incar.hiv.multiplier
        else:
            hiv_multiplier = 1

        # agent is incarcerated
        if self.active:
            self.time -= 1

            # FREE AGENT
            if self.time == 0:
                self.add_new_release_agent(self.agent)
                self.active = False

                # become high risk on release
                if (
                    not self.agent.high_risk.active and model.params.features.high_risk
                ):  # If behavioral treatment on and agent HIV, ignore HR period.
                    self.agent.high_risk.become_high_risk()
                    for bond in self.agent.location.params.high_risk.partnership_types:
                        self.agent.mean_num_partners[
                            bond
                        ] += self.agent.location.params.high_risk.partner_scale
                        self.agent.target_partners[bond] = utils.poisson(
                            self.agent.mean_num_partners[bond], model.np_random
                        )
                    model.pop.update_partnerability(self.agent)

                # does agent stay on haart
                if hiv_bool:
                    if self.agent.haart.active:
                        if (
                            model.run_random.random()
                            <= self.agent.location.params.incar.haart.discontinue
                        ):  # 12% remain surpressed
                            self.agent.haart.active = False
                            self.agent.haart.adherence = 0

                        # END FORCE

        # should the agent become incarcerated?
        elif model.run_random.random() < (
            self.agent.location.params.demographics[self.agent.race][
                self.agent.sex_type
            ].incar.prob
            * hiv_multiplier
            * model.calibration.incarceration
        ):
            incar_duration = self.agent.location.params.demographics[self.agent.race][
                self.agent.sex_type
            ].incar.duration.prob

            bin = current_p_value = 1
            p = model.run_random.random()
            while p >= current_p_value:
                current_p_value += incar_duration[bin].prob
                bin += 1

            timestay = model.run_random.randint(
                incar_duration[bin].min, incar_duration[bin].max
            )

            self.active = True
            self.time = timestay

            if hiv_bool:
                if not self.agent.hiv_dx:
                    if (
                        model.run_random.random()
                        < self.agent.location.params.incar.hiv.dx
                    ):
                        self.agent.hiv_dx = True
                else:  # Then tested and HIV, check to enroll in ART
                    if (
                        model.run_random.random()
                        < self.agent.location.params.incar.haart.prob
                    ):
                        if (
                            model.run_random.random()
                            < self.agent.location.params.incar.haart.adherence
                        ):
                            adherence = 5
                        else:
                            adherence = model.run_random.randint(1, 4)

                        # Add agent to HAART class set, update agent params
                        self.agent.haart.active = True
                        self.agent.haart.adherence = adherence
                        self.agent.haart.time = model.time

            # PUT PARTNERS IN HIGH RISK
            for bond in self.agent.location.params.high_risk.partnership_types:
                for partner in self.agent.partners[bond]:
                    if not partner.high_risk.active and model.params.features.high_risk:
                        if (
                            model.run_random.random()
                            < partner.location.params.high_risk.prob
                        ):
                            partner.high_risk.become_high_risk()
