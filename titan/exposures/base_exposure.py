from typing import List, Dict

from .. import agent
from .. import population
from .. import model


class BaseExposure:

    name: str = ""
    """Name of exposure in the params file.  Also used to name the attribute in Agent"""

    stats: List[str] = []
    """List of names of stats that come from this exposure (e.g. hiv.dx)"""

    def __init__(self, agent: "agent.Agent"):
        self.active = False
        self.agent = agent

    @classmethod
    def init_class(cls, params):
        """
        Initialize any class level attributes (such as setting counters to zero). Called on every active exposure on population initialization.

        args:
            params: parameters for this population
        """
        pass

    def init_agent(self, pop: "population.Population", time: int):
        """
        Initialize the agent for this exposure during population initialization (`Population.create_agent`).  Called only on exposures that are enabled per the params.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        pass

    def update_agent(self, model: "model.TITAN"):
        """
        Update the agent for this exposure for a time step.  Called once per time step in `TITAN.update_all_agents`.

        args:
            model: the instance of TITAN currently being run
        """
        pass

    @classmethod
    def add_agent(cls, agent: "agent.Agent"):
        """
        Add an agent to the class (not instance).  This can be useful if tracking population level statistics or groups, such as counts or newly active agents.

        This method is not called from anywhere in the model, but creates a cohesive api with `remove_agent`, which is called from `Population.remove_agent`.

        args:
            agent: the agent to add to the class attributes
        """
        pass

    @classmethod
    def remove_agent(cls, agent: "agent.Agent"):
        """
        Remove an agent from the class (not instance).  This can be useful if tracking population level statistics or groups, such as counts.

        This method is called from `Population.remove_agent`, but may also need to be called within the feature if an agent transitions from `active == True` to `active == False`.

        args:
            agent: the agent to remove from the class attributes
        """
        pass

    def set_stats(self, stats: Dict[str, int], time: int):
        """
        Update the `stats` dictionary passed for this agent.  Called from `output.get_stats` for each enabled exposure in the model.

        The stats to be updated must be declared in the class attribute `stats` to make sure the dictionary has the expected keys/counter value initialized.

        args:
            stats: the dictionary to update with this agent's feature statistics
            time: the time step of the model when the stats are set
        """
        pass

    @staticmethod
    def expose(
        model: "model.TITAN",
        interaction: str,
        rel: "agent.Relationship",
        num_acts: int,
    ):
        """
        Expose a relationship to the exposure for a number of acts for a specific interaction type.  Typically, this determines if the exposure can cause conversion/change in one of the agents, then if so, determines the probability of that and then converts the succeptible agent.

        args:
            model: The running model
            interaction: The type of interaction (e.g. sex, injection)
            rel: The relationship where the interaction is occuring
            num_acts: The number of acts of that interaction
        """
        pass

    def get_transmission_probability(
        self,
        model: "model.TITAN",
        interaction: str,
        partner: "agent.Agent",
        num_acts: int,
    ) -> float:
        """
        Determines the probability of a transmission event from agent to partner based on
            interaction type.

        This is not called from anywhere else in the model, but is recommended as a way to structure the exposure method into getting the transmission probability, then doing the conversion.

        args:
            model: The running model
            interaction : The type of interaction (e.g. `sex`, `injection`)
            partner: The agent's partner
            num_acts: The number of acts where exposure occured

        returns:
            probability of transmission from agent to partner
        """
        return 0.0

    def convert(self, model: "model.TITAN"):
        """
        Convert the agent to the exposure (i.e. become active).

        args:
            model: The running model
        """
        pass

    def diagnose(self, model: "model.TITAN"):
        """
        Diagnose the agent with the exposure (if applicable).

        args:
            model: The running model
        """
        pass
