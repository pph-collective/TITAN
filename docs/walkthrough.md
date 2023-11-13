# Walkthrough
For documentation on `run_titan.py`, as well as information for running the model interactively, see the [Running TITAN](running.md).

## Parameter setup
TITAN uses the [`paraml`](https://pypi.org/project/paraml/) package for parameter management using YAML. In TITAN, parameters are sub-divided, with each parameter set referring to a main model parameter type (e.g. [classes](https://pph-collective.github.io/titan-params-app/#/params#classes-1) or agent [demographics](https://pph-collective.github.io/titan-params-app/#/params#demographics-1)), a model feature (e.g. [syringe services](https://pph-collective.github.io/titan-params-app/#/params#syringe_services-1)), or an exposure (e.g. [hiv](https://pph-collective.github.io/titan-params-app/#/params#hiv-1)). Parameters may be set up in one YAML file or a folder of YAML files.

This guide will walk through creating a params file to augment the existing Atlanta setting, creating a small population to practice running. In order to best follow the next steps of the guide, users are encouraged to set this up as a single parameter file entitled `my_params.yml` in your working directory.

### Model parameters
To start, we want to edit the parameter file to run a smaller model over less time so that it will run quickly locally. To do this, we will update [model parameters](https://pph-collective.github.io/titan-params-app/#/params#model-1) to change population size and number of timesteps. In addition, we can set the random seed to provide reproducibility. To achieve this, add the following to your `my_params.yml` file (without comments):
```
model:
  seed: # random seed
    ppl: 1234  # set random seed for population initiation
    run: 1234  # set random seed for the run
  num_pop: 100  # create small population
  time:
    num_steps: 12  # number of steps in model run
    burn_steps: 0  # number of "burn-in" steps (negative time)
```

### Demographic parameters
Your parameter file can also be used to update [demographic](https://pph-collective.github.io/titan-params-app/#/params#demographics-1) information and distributions for agents in the model. Here, we will change the percentage of agents of each race (`white` and `black`) in the model. Add the following to `my_params.yml` (or play around with the numbers yourself! Just make sure they add up to 1):
```
demographics:
  black:
    ppl: 0.5
  white:
    ppl: 0.5
```

### Adding a random trial
Features can also be added to the model via params. To turn on the random trial feature, first we need to activate it:
```
features:
  random_trial: true
```
We can then change the parameters for the trial. Let's make it select agents randomly:
```
random_trial:
  choice: all
```

### Other parameters
Full reference for all parameters in TITAN can be found in our [parameter web app](https://pph-collective.github.io/titan-params-app).

## Running your params file
To run your new params file, you can simply use the command:
```
run_titan -p my_params.yml -S atlanta
```
This will save the results of your model in a newly-made `results` directory, or overwrite previous results in an existing `results` directory.

## Sweeping parameters
TITAN can also sweep over a set of parameters, either stepwise from the command
line or from a CSV defining a set of parameters. For either method, use 
`.`-separated parameters. For example, you might change the probability of
components being treated in the random trial module using [random_trial.prob](https://pph-collective.github.io/titan-params-app/#/params#random_trial-1).

To run this from the command line, use `param:start:stop:step`. For example, to
vary the random trial probability from 0.5 to 1 with steps of size 0.1, you can
run:
```
run_titan -p my_params.yml -S atlanta -w random_trial.prob:0.5:1.0:0.1
```

To sweep from values with a file, create a CSV with columns named as the `.`-separated parameters. For example, you might change the probability of componenets being treated in the random trial module by using the column name `random_trial.prob` in a `sweep_val.csv` file and populating rows with random numbers between 0 and 1. To run this sweep over the first 10 rows of your CSV, you would use the command:
```
run_titan -p my_params.yml -S atlanta -W sweep_file.csv -R 1:10
```
In either of the above examples, your results directory now contains a report with outputs using all 10 values, as well as a file linking each run's id with its value for the probability.