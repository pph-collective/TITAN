from .. import model
from .. import agent


class BaseInteraction:

    name: str = ""
    """Name of interaction in the params file."""

    @classmethod
    def interact(cls, model: "model.TITAN", rel: "agent.Relationship"):
        """
        Given a model and a relation, have the agents in the relationship interact for a time step.

        args:
            model: The running model
            rel: The relationship where interaction is happening
        """
        num_acts = cls.get_num_acts(model, rel)

        if num_acts < 1:
            return

        for exposure in model.exposures:
            if model.time >= model.params[exposure.name].start_time:
                exposure.expose(model, cls.name, rel, num_acts)

    @classmethod
    def get_num_acts(cls, model: "model.TITAN", rel: "agent.Relationship") -> int:
        return 0
