from typing import List, Dict, Optional, Set

from . import base_exposure
from .. import agent
from .. import population
from .. import model
from .. import utils


class HIV(base_exposure.BaseExposure):

    name: str = "hiv"
    stats: List[str] = ["hiv", "hiv_dx", "hiv_aids", "hiv_new", "hiv_dx_new"]
    """
        HIV collects the following stats:

        * hiv - number of agents with active hiv
        * hiv_dx - number of agents with diagnosed hiv
        * hiv_aids - number of agents with aids
        * hiv_new - number of agents converted to hiv this timestep
        * hiv_dx_new - number of agents with diagnosed with hiv this timestep
    """

    dx_counts: Dict[str, Dict[str, int]] = {}
    """Counts of diagnosed agents by race and sex_type"""

    agents: Set["agent.Agent"] = set()
    """Agents with active hiv"""

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

        Based on demographic params for the agent, stochastically determine if hiv is active, and if active, at what past time point was the agent converted, if the agent is diagnosed, and if the agent has aids.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        agent_params = (
            self.agent.location.params.demographics[self.agent.race]
            .sex_type[self.agent.sex_type]
            .drug_type[self.agent.drug_type]
        )

        # HIV
        if (
            pop.pop_random.random() < agent_params.hiv.init
            and time >= pop.params.hiv.start_time
        ):
            self.active = True

            # if HIV, when did the agent convert? Random sample
            self.time = utils.safe_random_int(
                time - self.agent.location.params.hiv.max_init_time,
                time,
                pop.pop_random,
            )

            if pop.pop_random.random() < agent_params.hiv.aids.init:
                self.aids = True

            if pop.pop_random.random() < agent_params.hiv.dx.init:
                self.dx = True
                # agent was diagnosed at a random time between conversion and now
                self.dx_time = utils.safe_random_int(self.time, time, pop.pop_random)

            # add agent to class
            self.add_agent(self.agent)

    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this exposure for a time step.  Called once per time step in `TITAN.update_all_agents`. Agent level updates are done after population level updates.   Called on only exposures that are enabled per the params.

        If the agent is hiv+ and the model time is past the hiv start_time, determine if the agent becomes diagnosed if not yet diagnosed, and if the agent has progressed to aids.

        args:
            model: the instance of TITAN currently being run
        """
        if self.active and model.time >= model.params.hiv.start_time:
            if not self.dx:
                test_prob = (
                    self.agent.location.params.demographics[self.agent.race]
                    .sex_type[self.agent.sex_type]
                    .drug_type[self.agent.drug_type]
                    .hiv.dx.prob
                )

                # Rescale based on calibration param
                test_prob *= model.calibration.test_frequency

                if model.run_random.random() < test_prob:
                    self.diagnose(model)

            self.progress_to_aids(model)

    @classmethod
    def add_agent(cls, agent: "agent.Agent"):
        """
        Add an agent to the class (not instance).  This can be useful if tracking population level statistics or groups, such as counts or newly active agents.

        Add the agent to the `agents` set and if the agent is diagnosed, updated the `dx_counts`

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

        Remove the agent from the `agents` set and decrement the `dx_counts` if the agent was diagnosed.

        args:
            agent: the agent to remove from the class attributes
        """
        cls.agents.remove(agent)

        if agent.hiv.dx:  # type: ignore[attr-defined]
            cls.dx_counts[agent.race][agent.sex_type] -= 1

    def set_stats(self, stats: Dict[str, int], time: int):
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
        model: "model.TITAN",
        interaction: str,
        rel: "agent.Relationship",
        num_acts: int,
    ):
        """
        Expose a relationship to the exposure for a number of acts of a specific interaction type.  Typically, this is determining if the exposure can cause conversion/change in one of the agents, then if so determining the probability of that and then converting the succeptible agent.

        For hiv, one agent must be active and the other not for an exposure to cause conversion.

        args:
            model: The running model
            interaction: The type of interaction (e.g. sex, injection)
            rel: The relationship where the interaction is occuring
            num_acts: The number of acts of that interaction
        """
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
            model, interaction, partner, num_acts
        )

        if model.run_random.random() < p:
            # if agent HIV+ partner becomes HIV+
            partner.hiv.convert(model)  # type: ignore[attr-defined]

    def get_transmission_probability(
        self,
        model: "model.TITAN",
        interaction: str,
        partner: "agent.Agent",
        num_acts: int,
    ) -> float:
        """
        Determines the probability of an hiv transmission event from agent to partner based on
            interaction type and numer of acts. For sex acts, transmission probability is a
            function of the acquisition probability of the HIV- agent's sex role
            and the HIV+ agent's haart adherence, acute status, and dx risk reduction

        args:
            model: The running model
            interaction : "injection" or "sex"
            partner: HIV- Agent
            num_acts: The number of exposure interactions the agents had this time step

        returns:
            probability of transmission from agent to partner
        """
        # if this isn't an interaction where hiv can transmit, return 0% prob
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

    def convert(self, model: "model.TITAN"):
        """
        Agent becomes HIV agent. Update all appropriate attributes, sets and dictionaries.

        args:
            model: The model being run
        """
        if not self.active:
            self.active = True
            self.time = model.time
            self.agent.vaccine.active = False  # type: ignore[attr-defined]
            self.add_agent(self.agent)

        if self.agent.prep.active:  # type: ignore[attr-defined]
            self.agent.prep.progress(model, force=True)  # type: ignore[attr-defined]

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
            hiv_duration = time - self.time

            if self.agent.location.params.hiv.acute.duration >= hiv_duration >= 0:
                return True

        return False

    def progress_to_aids(self, model: "model.TITAN"):
        """
        Model the progression of HIV agents to AIDS agents

        args:
             model: the running model
        """
        aids_prob = self.agent.location.params.hiv.aids.prob
        p = self.agent.haart.aids_scale()  # type: ignore[attr-defined]

        if model.run_random.random() < p * aids_prob:
            self.aids = True
