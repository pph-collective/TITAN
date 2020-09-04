from .base_feature import BaseFeature


class Prep(BaseFeature):

    # class level attributes to track all Prep agents
    counts = None
    new_agents = set()

    stats = ["numPrEP", "newNumPrEP", "injectable_prep", "oral_prep"]

    def __init__(self, agent):
        super().__init__(agent)
        # agent level attributes
        self.active = False
        self.adherence = 0
        self.type = ""
        self.load = 0.0
        self.last_dose = 0

        self.awareness = False
        self.opinion = 0.00

        # pca here or elsewhere?

    def update_agent(self, agent, model):
        if not agent.hiv:
            if model.time >= agent.location.params.prep.start_time:
                if agent.prep.active:
                    self.progress(agent, model)
                elif self.eligible(agent):
                    self.initiate(agent, model)

    @classmethod
    def update_pop(cls, model):
        # population is updated before agents, so clear set at the beginning of updates
        cls.new_agents = set()

    @classmethod
    def add_agent(cls, agent):
        # set up if this is the first time being called
        if cls.counts is None:
            cls.init_class(agent.location.params)

        cls.counts[agent.race] += 1
        cls.new_agents.add(agent)

    @classmethod
    def remove_agent(cls, agent):
        cls.counts[agent.race] -= 10

    def set_stats(self, stats, agent):
        if agent in self.new_agents:
            stats["newNumPrEP"] += 1

        if self.active:
            stats["numPrEP"] += 1
            if self.type == "Inj":
                stats["injectable_prep"] += 1
            elif self.type == "Oral":
                stats["oral_prep"] += 1

    ## HELPER FUNCTIONS

    @classmethod
    def init_class(cls, params):
        cls.counts = {race: 0 for race in params.classes.races}

    def initiate(self, agent, model, force: bool = False):
        """
        Place agents onto PrEP treatment. PrEP treatment assumes that the agent knows their HIV status is negative.

        args:
            agent : agent being updated
            force : whether to force the agent to enroll instead of using the appropriate algorithm per the prep params
        """
        # Prep only valid for agents not on prep and are HIV negative
        if self.active or agent.hiv:
            return

        if self.counts is None:
            self.init_class(agent.location.params)

        if force:
            self.enroll(agent, model.run_random)
        else:
            if "Racial" in agent.location.params.prep.target_model:
                num_prep_agents = self.counts[agent.race]
                all_hiv_agents = model.pop.hiv_agents.members
                all_race = {a for a in model.pop.all_agents if a.race == agent.race}

                num_hiv_agents = len(all_hiv_agents & all_race)
                target_prep = (
                    len(all_race) - num_hiv_agents
                ) * agent.location.params.demographcis[agent.race][
                    agent.sex_type
                ].prep.coverage
            else:
                num_prep_agents = sum(self.counts.values())
                target_prep = int(
                    (
                        model.pop.all_agents.num_members()
                        - model.pop.hiv_agents.num_members()
                    )
                    * agent.location.params.prep.target
                )

            if (
                num_prep_agents < target_prep
                and model.time >= agent.location.params.prep.start_time
                and self.eligible(agent)
            ):
                self.enroll(agent, model.run_random)

    def enroll(self, agent, rand_gen):
        """
        Enroll an agent in PrEP

        args:
            rand_gen: random number generator
        """
        params = agent.location.params

        self.active = True
        self.load = params.prep.peak_load
        self.last_dose = 0

        if (
            rand_gen.random()
            < params.demographics[agent.race][agent.sex_type].prep.adherence
        ):
            self.adherence = 1
        else:
            self.adherence = 0

        if "Inj" in params.prep.type and "Oral" in params.prep.type:
            if rand_gen.random() < params.prep.lai.prob:
                self.type = "Inj"
            else:
                self.type = "Oral"

        else:
            self.type = params.prep.type[0]

        self.add_agent(agent)

    def progress(self, agent, model, force: bool = False):
        """
        Update agent's PrEP status and discontinue stochastically or if `force` is True

        args:
            agent: agent being updated
            force: whether to force discontinuation of PrEP
        """
        if force:
            self.remove_agent(agent)
            return

        if (
            model.run_random.random()
            < agent.location.params.demographics[agent.race][
                agent.sex_type
            ].prep.discontinue
            and self.type == "Oral"
        ):
            self.discontinue(agent)

        if self.type == "Inj":
            self.update_load(agent)

    def update_load(self, agent):
        """
        Determine and update load of PrEP concentration in agent.
        """
        params = agent.location.params

        self.last_dose += 1
        if self.last_dose > params.model.time.steps_per_year:
            self.remove_agent(agent)
            self.load = 0.0
        else:
            annualized_last_dose = self.last_dose / params.model.time.steps_per_year
            annualized_half_life = params.prep.half_life / 365
            self.load = params.prep.peak_load * (
                (0.5) ** (annualized_last_dose / annualized_half_life)
            )

    def discontinue(self, agent):
        self.active = False
        self.type = ""
        self.load = 0.0
        self.last_dose = 0

        self.remove_agent(agent)

    def eligible(self, agent):
        """
        Determine if an agent is eligible for PrEP

        returns:
            prep
        """
        target_model = agent.location.params.prep.target_model
        gender = agent.location.params.classes.sex_types[agent.sex_type].gender

        if self.active or agent.vaccine.active:
            return False

        all_eligible_models = {"Allcomers", "Racial"}

        if all_eligible_models.intersection(target_model):
            return True

        if "cdc_women" in target_model:
            if gender == "F":
                if self.cdc_eligible(agent):
                    return True

        if "cdc_msm" in target_model:
            if gender == "M" and self.cdc_eligible(agent):
                return True

        if "pwid_sex" in target_model:
            if agent.drug_type == "Inj" and self.cdc_eligible(agent):
                return True

        if "pwid" in target_model:
            if agent.drug_type == "Inj":
                return True

        if "ssp_sex" in target_model:
            if agent.ssp and self.cdc_eligible(agent):
                return True

        if "ssp" in target_model:
            if agent.ssp:
                return True

        return False

    @staticmethod
    def cdc_eligible(agent) -> bool:
        """
        Determine agent eligibility for PrEP under CDC criteria

        returns:
            cdc eligibility
        """
        if self.is_msm():
            return True

        ongoing_duration = agent.location.params.partnership.ongoing_duration
        for rel in agent.relationships:
            partner = rel.get_partner(self)
            if rel.duration > ongoing_duration and partner.hiv_dx:
                return True

            if partner.drug_type == "Inj" or partner.is_msm():
                return True

        return False
