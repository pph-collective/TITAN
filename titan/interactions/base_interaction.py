from .. import model
from .. import agent


class BaseInteraction:

    name: str = ""
    """Name of interaction in the params file."""

    @staticmethod
    def interact(model: "model.HIVModel", rel: "agent.Relationship") -> bool:
        """
        Given a model and a relation, have the agents in the relationship interact for a time step.

        args:
            model: The running model
            rel: The relationship where interaction is happening

        returns:
            whether the agents interacted
        """
        pass
