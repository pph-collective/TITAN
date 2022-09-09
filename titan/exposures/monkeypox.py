from typing import List, Dict, Optional, Set

from . import base_exposure
from .. import agent
from .. import population
from .. import model
from .. import utils


class MonkeyPox(base_exposure.BaseExposure):

    name: str = "monkeypox"
    stats: List[str] = [
        "monkeypox",
        "monkeypox_dx",
        "monkeypox_new",
        "monkeypox_dx_new",
    ]
    """
        MonkeyPox collects the following stats:

        * monkeypox - number of agents with hx of monkeypox
        * monkeypox_dx - number of agents ever diagnosed with monkeypox
        * monkeypox_new - number of agents converted to monkeypox this timestep
        * monkeypox_dx_new - number of agents with diagnosed with monkeypox this timestep
    """

    dx_counts: Dict[str, Dict[str, int]] = {}
    """Counts of diagnosed agents by race and sex_type"""

    agents: Set["agent.Agent"] = set()
    """Agents who have ever had monkeypox"""

    def __init__(self, agent: "agent.Agent"):
        super().__init__(agent)

        self.active = False
        self.time: Optional[int] = None
        self.dx = False
        self.dx_time: Optional[int] = None

    @classmethod
    def init_class(cls, params):
        """
        Initialize any diagnosis counts and the agents set.

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
        Initialize the agent for this exposure during population initialization (`Population.create_agent`).  Called on only exposures that are enabled per the params.

        Based on demographic params for the agent, stochastically determine if monkeypox is active, and if active, at what past time point was the agent converted, if the agent is diagnosed, and if the agent has aids.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        agent_params = (
            self.agent.location.params.demographics[self.agent.race]
            .sex_type[self.agent.sex_type]
            .drug_type[self.agent.drug_type]
        )

        # Monkeypox
        if (
            pop.pop_random.random() < agent_params.monkeypox.init
            and time >= pop.params.monkeypox.start_time
        ):
            self.active = True

            # if monkeypox, when did the agent convert? Random sample
            self.time = utils.safe_random_int(
                time - self.agent.location.params.monkeypox.max_init_time,
                time,
                pop.pop_random,
            )

            if pop.pop_random.random() < agent_params.monkeypox.dx.init:
                self.dx = True
                # agent was diagnosed at a random time between conversion and now
                self.dx_time = utils.safe_random_int(self.time, time, pop.pop_random)

            # add agent to class
            self.add_agent(self.agent)

    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this exposure for a time step.  Called once per time step in `TITAN.update_all_agents`. Agent level updates are done after population level updates. Called on only exposures that are enabled per the params.

        If the agent is monkeypox+ and the model time is past the monkeypox start_time, determine if the agent becomes diagnosed if not yet diagnosed.

        args:
            model: the instance of TITAN currently being run
        """
        if self.active and model.time >= model.params.monkeypox.start_time:
            if not self.dx:
                test_prob = (
                    self.agent.location.params.demographics[self.agent.race]
                    .sex_type[self.agent.sex_type]
                    .drug_type[self.agent.drug_type]
                    .monkeypox.dx.prob
                )

                # Rescale based on calibration param
                test_prob *= model.calibration.test_frequency

                if model.run_random.random() < test_prob:
                    self.diagnose(model)

    @classmethod
    def add_agent(cls, agent: "agent.Agent"):
        """
        Add an agent to the class (not instance).  This can be useful if tracking population level statistics or groups, such as counts or newly active agents.

        Add the agent to the `agents` set and if the agent is diagnosed, updated the `dx_counts`

        args:
            agent: the agent to add to the class attributes
        """
        cls.agents.add(agent)

        if agent.monkeypox.dx:  # type: ignore[attr-defined]
            cls.dx_counts[agent.race][agent.sex_type] += 1

    @classmethod
    def remove_agent(cls, agent: "agent.Agent"):
        """
        Remove an agent from the class (not instance).  This can be useful if tracking population level statistics or groups, such as counts.

        Remove the agent from the `agents` set and decrement the `dx_counts` if the agent was diagnosed.

        args:
            agent: the agent to remove from the class attributes
        """
        cls.agents.remove(agent)

        if agent.monkeypox.dx:  # type: ignore[attr-defined]
            cls.dx_counts[agent.race][agent.sex_type] -= 1

    def set_stats(self, stats: Dict[str, int], time: int):
        if self.active:
            stats["monkeypox"] += 1
            if self.time == time:
                stats["monkeypox_new"] += 1
            if self.dx:
                stats["monkeypox_dx"] += 1
                if self.dx_time == time:
                    stats["monkeypox_dx_new"] += 1

    @staticmethod
    def expose(
        model: "model.TITAN",
        interaction: str,
        rel: "agent.Relationship",
        num_acts: int,
    ):
        """
        Expose a relationship to the exposure for a number of acts of a specific interaction type.  Typically, this is determining if the exposure can cause conversion/change in one of the agents, then if so determining the probability of that and then converting the succeptible agent.

        For monkeypox, one agent must be active and the other not for an exposure to cause conversion.

        args:
            model: The running model
            interaction: The type of interaction (e.g. sex, injection)
            rel: The relationship where the interaction is occuring
            num_acts: The number of acts of that interaction
        """
        # Agent 1 is monkeypox+, Agent 2 is not, Agent 2 is succept
        if rel.agent1.monkeypox.active and not rel.agent2.monkeypox.active:  # type: ignore[attr-defined]
            agent = rel.agent1
            partner = rel.agent2
        # If Agent 2 is monkeypox+ and Agent 1 is not, Agent 1 is succept
        elif not rel.agent1.monkeypox.active and rel.agent2.monkeypox.active:  # type: ignore[attr-defined]
            agent = rel.agent2
            partner = rel.agent1
        else:  # neither agent is monkeypox+ or both are
            return

        p = agent.monkeypox.get_transmission_probability(  # type: ignore[attr-defined]
            model, interaction, partner, num_acts
        )

        if model.run_random.random() < p:
            # if agent monkeypox+ partner becomes monkeypox+
            partner.monkeypox.convert(model)  # type: ignore[attr-defined]

    def get_transmission_probability(
        self,
        model: "model.TITAN",
        interaction: str,
        partner: "agent.Agent",
        num_acts: int,
    ) -> float:
        """
        Determines the probability of an monkeypox transmission event from agent to partner based on
            interaction type and numer of acts. For sex acts, transmission probability is a
            function of the acquisition probability of the monkeypox- agent's sex role
            and the monkeypox+ agent's haart adherence, acute status, and dx risk reduction

        args:
            model: The running model
            interaction : "injection" or "sex"
            partner: monkeypox- Agent
            num_acts: The number of exposure interactions the agents had this time step

        returns:
            probability of transmission from agent to partner
        """
        # if this isn't an interaction where monkeypox can transmit, return 0% prob
        if interaction not in ("sex") or not self.get_acute_status(model.time):
            return 0.0

        # Logic for if needle or sex type interaction
        p: float

        # get partner's sex role during acts
        partner_sex_role = "versatile"

        # get probability of sex acquisition given monkeypox- partner's position
        p = partner.location.params.partnership.sex.acquisition[partner.sex_type][
            partner_sex_role
        ]

        # feature specific risk adjustment
        for feature in model.features:
            agent_feature = getattr(self.agent, feature.name)
            p *= agent_feature.get_transmission_risk_multiplier(self.time, interaction)

            partner_feature = getattr(partner, feature.name)
            p *= partner_feature.get_acquisition_risk_multiplier(self.time, interaction)

        # Scaling parameter for positively identified monkeypox agents
        if self.dx:
            p *= 1 - self.agent.location.params.monkeypox.dx.risk_reduction[interaction]

        # Racial calibration parameter to attain proper race incidence disparity
        p *= partner.location.params.demographics[partner.race].monkeypox.transmission

        # Scaling parameter for per act transmission.
        p *= model.calibration.acquisition

        return utils.total_probability(p, num_acts)

    def convert(self, model: "model.TITAN"):
        """
        Agent becomes monkeypox agent. Update all appropriate attributes, sets and dictionaries.

        args:
            model: The model being run
        """
        if not self.active:
            self.active = True
            self.time = model.time
            self.agent.vaccine.active = False  # type: ignore[attr-defined]
            self.add_agent(self.agent)

    def diagnose(self, model: "model.TITAN"):
        """
        Mark the agent as diagnosed.

        args:
             model: the running model
        """
        self.dx = True
        self.dx_time = model.time
        self.add_agent(self.agent)

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
            monkeypox_duration = time - self.time

            if (
                self.agent.location.params.monkeypox.acute.duration
                >= monkeypox_duration
                >= 0
            ):
                return True

        return False
