from typing import Dict

from .base_feature import BaseFeature


class HAART(BaseFeature):
    """
    Highly Active Antiretroviral Theray (HAART) is a treatment regimen.
    """

    name = "haart"
    stats = ["haart"]
    """
        HAART collects the following stats:

        * haart - number of agents with active haart
    """

    counts = None

    def __init__(self, agent: "Agent"):
        super().__init__(agent)

        self.active = False
        self.time = 0
        self.adherence = 0

    def init_agent(self, pop: "Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        An Agent can only be initialized with HAART if they are HIV+ and diagnosed.  They are randomly assigned to HAART with a probability based on demographics and then if assigned to haart, assigned an adherence based on those same demographic params.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
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

    def update_agent(self, model: "HIVModel"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `HIVModel.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        Account for HIV treatment through highly active antiretroviral therapy (HAART). HAART was implemented in 1996, hence, there is treatment only after 1996. HIV treatment assumes that the agent knows their HIV+ status (`dx` is True).

        args:
            model: the instance of HIVModel currently being run
        """
        if self.agent.hiv and model.time >= model.params.hiv.start_time:
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
                        num_haart_agents = self.counts[self.agent.race][
                            self.agent.sex_type
                        ]

                        if num_haart_agents < (haart_params.prob * num_dx_agents):
                            self.initiate(model)
                    else:
                        if model.run_random.random() < (
                            haart_params.prob * model.calibration.haart.coverage
                        ):
                            self.initiate(model)
                # Go off HAART
                elif (
                    self.active and model.run_random.random() < haart_params.discontinue
                ):
                    self.active = False
                    self.adherence = 0
                    self.time = 0
                    self.remove_agent(self.agent)

    @classmethod
    def add_agent(cls, agent: "Agent"):
        """
        Add an agent to the class (not instance).

        Increments `counts` or haart agents by race and sex_type for the given agent.

        args:
            agent: the agent to add to the class attributes
        """
        # set up if this is the first time being called
        if cls.counts is None:
            cls.init_class(agent.location.params)

        cls.counts[agent.race][agent.sex_type] += 1

    def remove_agent(cls, agent: "Agent"):
        """
        Remove an agent from the class (not instance).

        Decrements `counts` or haart agents by race and sex_type for the given agent.

        args:
            agent: the agent to remove from the class attributes
        """
        cls.counts[agent.race][agent.sex_type] -= 1

    def set_stats(self, stats: Dict[str, int]):
        if self.active:
            stats["haart"] += 1

    # =========== HELPER METHODS ============

    @classmethod
    def init_class(cls, params: "DotMap"):
        """
        Initialize the counts dictionary for the races and sex_types in the model.

        args:
            params: the model params
        """
        cls.counts = {
            race: {sex_type: 0 for sex_type in params.classes.sex_types}
            for race in params.classes.races
        }

    def initiate(self, model: "HIVModel"):
        """
        Initiate an agent with HAART and add them to the population.

        args:
            model: the instance of HIVModel currently being run
        """
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
