## Agent Interactions

Interactions from [`params.classes.bond_types.acts_allowed`](https://pph-collective.github.io/titan-params-app/#/params#classes-1) are implemented in standalone files that implement the interface from `interactions.BaseInteraction`.  This allows for the logic related to an interaction type to be consolidated into one place and make incorporating a new interaction type as simple as possible.

To add a new interaction:

* Add the param to `param.classes.bond_types.acts_allowed`
* Add a file to `interactions/` which creates a class that is a sub-class of `BaseInteraction`
    * Implement the methods of `BaseInteraction` which are needed for this interaction
    * Not all methods are needed for all interactions (see below for details on the methods)
* Re-export the feature from `interactions/__init__.py`
* Add tests in `tests/interactions/`
* Add it to the docs in `docs/api/interactions/` and to the nav in `mkdocs.yml`

The `TITAN` class uses sub-classes of `BaseInteraction` to initialize the object/call methods as appropriate.

::: titan.interactions.base_interaction
