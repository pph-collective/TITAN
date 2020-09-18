class BaseFeature:
    """
    Interface class for an agent-oriented feature of the TITAN model.  It is intended to be
    inherited by feature classes to ensure that the expected methods/fields exist.

    The class takes advantage of both instance and class methods and attributes. The instance
    methods/attributes are used at the agent level, whereas the class methods/attributes are
    used at the population level.

    A class method can be called on an instance of an object, but it doesn't have access
    to that instance (e.g. agent.feature.add_agent has to be passed agent because it is a class
    method).
    """

    """Name of feature in the params file.  Also used to name the attribute in Agent"""
    name = ""

    """List of names of stats that come from this feature (e.g. numFeat)"""
    stats = []

    def __init__(self, agent):
        """
        Constructor for an instance of the feature.  This is called from within `Agent.__init__` and passes the agent to the feature to create a two way binding.  All features must have the attributes of `active` and `agent`.  By default `active` is false and `agent` is the passed agent.

        Called on all features whether or not param is enabled to make sure that references to other features within a feature do not cause errors (e.g. incar referring to prep, but prep isn't on).

        args:
            agent: The agent this feature instance is attached to.
        """
        self.active = False
        self.agent = agent

    def update_agent(self, model):
        """
        Update the agent for this feature for a time step.  Called once per time step in `HIVModel.update_all_agents`. Agent level updates are done after population level updates.   Called on only features that are enabled per the params.

        args:
            model: the instance of HIVModel currently being run
        """
        pass

    def init_agent(self, pop, time):
        """
        Initialize the agent for this feature during population initialization (`Population.create_agent`).  Called on only features that are enabled per the params.

        args:
            pop: the population this agent is a part of
            time: the current time step
        """
        pass

    @classmethod
    def add_agent(cls, agent):
        """
        Add an agent to the class (not instance).  This can be useful if tracking population level statistics or groups, such as counts or newly active agents.

        This method is not called from anywhere in the model, but creates a cohesive api with `remove_agent`, which is called from `Population.remove_agent`.

        args:
            agent: the agent to add to the class attributes
        """
        pass

    @classmethod
    def remove_agent(cls, agent):
        """
        Remove an agent from the class (not instance).  This can be useful if tracking population level statistics or groups, such as counts.

        This method is called from `Population.remove_agent`, but may also need to be called within the feature if an agent transitions from `active == True` to `active == False`.

        args:
            agent: the agent to remove from the class attributes
        """
        pass

    @classmethod
    def update_pop(cls, model):
        """
        Update the feature for the entire population (class method).  This is useful for initializing class level trackers that need to be reset each time step, or if enabling a feature for agents needs to be evaluated within the context of the full population (limited slots, or similar).

        This is called in `HIVModel.update_all_agents` before agent-level updates are made.

        args:
            model: the instance of HIVModel currently being run
        """
        pass

    def set_stats(self, stats):
        """
        Update the `stats` dictionary passed for this agent.  Called from `output.get_stats` for each enabled feature in the model.

        The stats to be updated must be declared in the class attribute `stats` to make sure the dictionary has the expected keys/counter value initialized.

        args:
            stats: the dictionary to update with this agent's feature statistics
        """
        pass
