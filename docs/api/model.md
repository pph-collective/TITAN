## TITAN

The `TITAN` class is used to model agent interactions as they progress through time.  The model can be run on an existing `Population` or it can create a `Population` during its construction.  The most common entry point to the model is `run`, which will run the model for all time steps.  To run step by step, `step` can be iterated through instead, just be sure to `reset_trackers` between `step`s.

::: titan.model.TITAN
