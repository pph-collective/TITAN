from typing import List, Dict, Optional, Set

from . import base_exposure
from .. import agent
from .. import population
from .. import model
from .. import probabilities as prob
from .. import utils


class HIV(base_exposure.BaseExposure):

    name: str = "hiv"
    """Name of exposure in the params file.  Also used to name the attribute in Agent"""

    stats: List[str] = ["hiv", "hiv_dx", "hiv_aids", "hiv_new", "hiv_dx_new"]
    """List of names of stats that come from this exposure (e.g. hiv.dx)"""

    dx_counts: Dict[str, Dict[str, int]] = {}
    agents: Set["agent.Agent"] = set()

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)

        self.active = False
        self.time: Optional[int] = None
        self.dx = False
        self.dx_time: Optional[int] = None
        self.aids = False

    @classmethod
    def init_class(cls, params):
        """
        Initialize any class level attributes (such as setting counters to zero). Called on every active feature on population initialization.

        args:
            params: parameters for this population
        """
        cls.dx_counts = {
            race: {so: 0 for so in params.classes.sex_types}
            for race in params.classes.races
        }
        cls.agents = set()

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        agent_params = self.agent.location.params.demographics[self.agent.race][
            self.agent.population
        ]

        # HIV
        if (
            pop.pop_random.random() < agent_params.hiv.init
            and time >= pop.params.hiv.init
        ):
            self.active = True
            self.add_agent(self.agent)

            # if HIV, when did the agent convert? Random sample
            self.time = pop.pop_random.randint(
                time - self.agent.location.params.hiv.max_init_time, time
            )

            if pop.pop_random.random() < agent_params.aids.init:
                self.aids = True

            if pop.pop_random.random() < agent_params.hiv.dx.init:
                self.dx = True
                # agent was diagnosed at a random time between conversion and now
                self.dx_time = pop.pop_random.randint(self.time, time)

    def update_agent(self, model: "model.HIVModel"):
        """
        Update the agent for this feature for a time step.  Called once per time step in `HIVModel.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        args:
            model: the instance of HIVModel currently being run
        """
        if self.active and model.time >= model.params.hiv.start_time:
            partner_tracing = self.agent.location.params.partner_tracing

            if not self.dx:
                test_prob = self.agent.location.params.demographics[self.agent.race][
                    self.agent.sex_type
                ].hiv.dx.prob

                # Rescale based on calibration param
                test_prob *= model.calibration.test_frequency

                if model.run_random.random() < test_prob:
                    self.diagnose(model)
                elif (
                    self.agent.partner_traced
                    and model.run_random.random() < partner_tracing.hiv.dx
                    and model.time > self.agent.trace_time
                ):
                    self.diagnose(model)

            if model.time >= self.agent.trace_time + partner_tracing.trace_duration:
                # agents can only be traced during a specified period after their partner is
                # diagnosed. If past this time, remove ability to trace.
                self.agent.partner_traced = False

            self.progress_to_aids(model)

    @classmethod
    def add_agent(cls, agent: "agent.Agent"):
        """
        Add an agent to the class (not instance).  This can be useful if tracking population level statistics or groups, such as counts or newly active agents.

        This method is not called from anywhere in the model, but creates a cohesive api with `remove_agent`, which is called from `Population.remove_agent`.

        args:
            agent: the agent to add to the class attributes
        """
        cls.agents.add(agent)

        if agent.hiv.dx:  # type: ignore[attr-defined]
            cls.dx_counts[agent.race][agent.sex_type] += 1

    @classmethod
    def remove_agent(cls, agent: "agent.Agent"):
        """
        Remove an agent from the class (not instance).  This can be useful if tracking population level statistics or groups, such as counts.

        This method is called from `Population.remove_agent`, but may also need to be called within the feature if an agent transitions from `active == True` to `active == False`.

        args:
            agent: the agent to remove from the class attributes
        """
        cls.agents.remove(agent)

        if agent.hiv.dx:  # type: ignore[attr-defined]
            cls.dx_counts[agent.race][agent.sex_type] -= 1

    def set_stats(self, stats: Dict[str, int], time: int):
        """
        Update the `stats` dictionary passed for this agent.  Called from `output.get_stats` for each enabled feature in the model.

        The stats to be updated must be declared in the class attribute `stats` to make sure the dictionary has the expected keys/counter value initialized.

        args:
            stats: the dictionary to update with this agent's feature statistics
            time: the time step of the model when the stats are set
        """
        if self.active:
            stats["hiv"] += 1
            if self.time == time:
                stats["hiv_new"] += 1
            if self.aids:
                stats["hiv_aids"] += 1
            if self.dx:
                stats["hiv_dx"] += 1
                if self.dx_time == time:
                    stats["hiv_dx_new"] += 1

    @staticmethod
    def expose(
        model: "model.HIVModel",
        interaction: str,
        rel: "agent.Relationship",
        num_acts: int,
    ):
        # Agent 1 is HIV+, Agent 2 is not, Agent 2 is succept
        if rel.agent1.hiv.active and not rel.agent2.hiv.active:  # type: ignore[attr-defined]
            agent = rel.agent1
            partner = rel.agent2
        # If Agent 2 is HIV and Agent 1 is not, Agent 1 is succept
        elif not rel.agent1.hiv.active and rel.agent2.hiv.active:  # type: ignore[attr-defined]
            agent = rel.agent2
            partner = rel.agent1
        else:  # neither agent is HIV or both are
            return

        p = agent.hiv.get_transmission_probability(  # type: ignore[attr-defined]
            model, "sex", partner, num_acts
        )

        if model.run_random.random() < p:
            # if agent HIV+ partner becomes HIV+
            partner.hiv.convert(model)  # type: ignore[attr-defined]

    def get_transmission_probability(
        self,
        model: "model.HIVModel",
        interaction: str,
        partner: "agent.Agent",
        num_acts: int,
    ) -> float:
        """
        Determines the probability of an hiv transmission event from agent to partner based on
            interaction type. For sex acts, transmission probability is a
            function of the acquisition probability of the HIV- agent's sex role
            and the HIV+ agent's haart adherence, acute status, and dx risk reduction

        args:
            interaction : "injection" or "sex"
            agent: HIV+ Agent
            partner: HIV- Agent

        returns:
            probability of transmission from agent to partner
        """
        # if this isn't an interaction where hiv can transmit, return 0.0
        if interaction not in ("injection", "sex"):
            return 0.0

        # Logic for if needle or sex type interaction
        p: float

        # get baseline probabilities
        if interaction == "injection":
            p = model.params.partnership.injection.transmission.base
        elif interaction == "sex":
            agent_sex_role = self.agent.sex_role
            partner_sex_role = partner.sex_role

            # get partner's sex role during acts
            if partner_sex_role == "versatile":  # versatile partner takes
                # "opposite" position of agent
                if agent_sex_role == "insertive":
                    partner_sex_role = "receptive"
                elif agent_sex_role == "receptive":
                    partner_sex_role = "insertive"
                else:
                    partner_sex_role = "versatile"  # if both versatile, can switch
                    # between receptive and insertive by act

            # get probability of sex acquisition given HIV- partner's position
            p = partner.location.params.partnership.sex.acquisition[partner.sex_type][
                partner_sex_role
            ]

        # feature specific risk adjustment
        for feature in model.features:
            agent_feature = getattr(self.agent, feature.name)
            p *= agent_feature.get_transmission_risk_multiplier(self.time, interaction)

            partner_feature = getattr(partner, feature.name)
            p *= partner_feature.get_acquisition_risk_multiplier(self.time, interaction)

        # Scaling parameter for acute HIV infections
        if self.get_acute_status(model.time):
            p *= self.agent.location.params.hiv.acute.infectivity

        # Scaling parameter for positively identified HIV agents
        if self.dx:
            p *= 1 - self.agent.location.params.hiv.dx.risk_reduction[interaction]

        # Racial calibration parameter to attain proper race incidence disparity
        p *= partner.location.params.demographics[partner.race].hiv.transmission

        # Scaling parameter for per act transmission.
        p *= model.calibration.acquisition

        return utils.total_probability(p, num_acts)

    def convert(self, model: "model.HIVModel"):
        """
        Agent becomes HIV agent. Update all appropriate list and dictionaries.

        args:
            agent: The agent being converted
        """
        if not self.active:
            self.active = True
            self.time = model.time
            self.agent.vaccine.active = False  # type: ignore[attr-defined]
            self.add_agent(self.agent)

        if self.agent.prep.active:  # type: ignore[attr-defined]
            self.agent.prep.progress(model, force=True)  # type: ignore[attr-defined]

    # ============================ HELPER METHODS ==============================

    def get_acute_status(self, time: int) -> bool:
        """
        Get acute status of agent at time

        args:
            time: The current time step

        returns:
            whether an agent is acute
        """
        if self.active and self.time is not None:
            hiv_duration = time - self.time

            if self.agent.location.params.hiv.acute.duration >= hiv_duration >= 0:
                return True

        return False

    def diagnose(self, model: "model.HIVModel"):
        """
        Mark the agent as diagnosed and trace their partners (if partner tracing enabled).

        args:
             model: the running model
        """
        # agent's location's params used throughout as that is the agent who
        # would be interacting with the service
        partner_tracing = self.agent.location.params.partner_tracing

        self.dx = True
        self.dx_time = model.time
        self.add_agent(self.agent)

        # do partner tracing if enabled
        if (
            model.params.features.partner_tracing
            and partner_tracing.start_time <= self.time < partner_tracing.stop_time
        ):
            # Determine if each partner is found via partner tracing
            for ptnr in self.agent.get_partners(partner_tracing.bond_type):
                if not ptnr.hiv.dx and model.run_random.random() < partner_tracing.prob:  # type: ignore[attr-defined]
                    ptnr.partner_traced = True
                    ptnr.trace_time = model.time

    def progress_to_aids(self, model: "model.HIVModel"):
        """
        Model the progression of HIV agents to AIDS agents

        args:
             model: the running model
        """
        p = prob.adherence_prob(self.agent.haart.adherence) if self.agent.haart.active else 1  # type: ignore[attr-defined]

        if model.run_random.random() < p * self.agent.location.params.hiv.aids.prob:
            self.aids = True
