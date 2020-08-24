## Population

The `Population` class is used to represent the population of agents the model is running on.  On construction, it stochastically creates the population described in the `params`.  At its core, it is a graph with nodes (`all_agents`) and edges (`relationships`), it can be formally backed by a NetworkX graph by enabling the graph in the prams file.  This allows for some graph-specific logic to be applied throughout the running of the model (e.g. trimming components, writing network statistics).

::: titan.population.Population

## Population Reading & Writing

!!! info "Released in v1.1.0"

Populations can be saved to file so that they can be analysed in detail or re-used in a future run.  `run_titan.py` allows this using the `--savepop [path]` option to save the population to the path, and the `--poppath [path]` option loads the population at the path.  The population is saved after creation, but before the model has run.

### Saving the Population

The population is represented as a series of csv files that save the attributes for the core entities (agents, relationships at this time).  The population can be saved with only `core` attributes (e.g. race, sex_type, hiv) or with `intervention` attributes (e.g. prep, vaccinated) as well.  `intervention` attributes are less likely to work as intended across versions of the model.

::: titan.population_io.write
    rendering:
        show_root_heading: true
        show_root_toc_entry: false
        heading_level: 4

### Reading in/using a Saved Population

A population can be created from the files saved, however, there is no validation done to ensure the params used when running on this population make sense/match what was originally used when creating it.  Some things may validly change (e.g. interventions, reports needed, seeds), but others may result in strange behavior if changed (e.g. race distribution, what classes are in use).

::: titan.population_io.read
    rendering:
        show_root_heading: true
        show_root_toc_entry: false
        heading_level: 4
