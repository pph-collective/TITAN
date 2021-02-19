## Running the Model

The model has a wrapper script called `run_titan.py` that makes running a full simulation easy.  TITAN can also be run from an interactive repl or a custom script.

### run_titan.py

To run the model, execute the `run_titan.py` program within the `/titan/` directory. See [TITAN params](https://pph-collective.github.io/titan-params-app) for documentation on how to set and use parameters.

Results of the model are generated and aggregated into the `/results/` directory by default. If the model is re-run, the existing results will be overwritten.

#### Usage

Below are the results of `python run_titan.py --help`.  It highlights all of the command line arguments that can be passed to the script.

```
usage: run_titan.py [-h] [-n [NMC]] [-S SETTING] -p PARAMS [-o OUTDIR]
                    [-b BASE] [-e] [--savepop] [--poppath POPPATH]
                    [-w SWEEP [SWEEP ...]] [-W SWEEPFILE] [-r ROWS] [-F]


Run TITAN model

optional arguments:
  -h, --help            show this help message and exit
  -n [NMC], --nMC [NMC]
                        number of monte carlo runs to complete
  -S SETTING, --setting SETTING
                        setting directory to use
  -p PARAMS, --params PARAMS
                        directory or file with params yaml(s)
  -o OUTDIR, --outdir OUTDIR
                        directory name to save results to
  -b BASE, --base BASE  whether to use base setting
  -e, --error           Error on unused parameters instead of warning
  --savepop             Save population after creation, but before model run.
  --poppath POPPATH     Path to saved population (directory or .tar.gz file)
  -w SWEEP [SWEEP ...], --sweep SWEEP [SWEEP ...]
                        Optional and repeatable definitions of numeric params
                        to sweep. Expected format is param:start:stop[:step]
  -W SWEEPFILE, --sweepfile SWEEPFILE
                        Optional. CSV file with param sweep definitions.
                        Header row must contain param paths, with data rows
                        containing values. If this is passed, any `-w` args
                        will be ignored.
  -r ROWS, --rows ROWS  Optional. Which data rows of sweepfile to use in
                        format start:stop.
  -F, --force           Run model even if number of sweeps exceeds 100
```

::: titan.run_titan.main
    rendering:
        show_root_heading: true
        show_root_toc_entry: false
        heading_level: 4


### Running Interactively

The model can also be run interactively in the repl.  Start a `python` session from the root directory of `TITAN`, and follow along!

We'll use the sample params file `tests/params/basic.yml` in all of these examples, but feel free to use a different one.

Here is how to perform the basic steps of running the model:
```python
from titan.parse_params import create_params
from titan.model import TITAN

outdir = 'results'

params = create_params(None, 'tests/params/basic.yml', outdir)
model = TITAN(params)
model.run(outdir)
```
This creates a params object using no setting (the `None`), our test params, and tells `create_params` to put our computed params file in a directory called `results`.

!!! note "the 'results' directory must already be created"

We then use those params to create our model, and run it.  We also have the model results saved to our `results` directory.

We should now see a `params.yml` in our 'results' directory, and some reports showing what happened at different timesteps in the model.

If we wanted to debug something, or look at a very specific metric that wasn't in our reports, we could instead step through the model one time-step at a time.

Resuming from our code above, here's how we could do that.
```python
model2 = TITAN(params)
start_time = 0
end_time = 10
for i in range(start_time, end_time):
    model2.time = i # update the model's time
    model2.step(outdir)

    # do some introspection here, like...
    print(model2.pop.haart_counts)

    # make sure the model state is reset for a new time step
    model2.reset_trackets()

```

If we want to write and read in a population instead of letting the model create one...

```python
from titan import population_io as pio
from titan.population import Population
from copy import deepcopy

# let's make a copy of our params and tinker with the population a bit
params2 = deepcopy(params)
params2.demographics.white.MSM.hiv.prob = 0.4

pop = Population(params2)
poppath = pio.write(pop, outdir)

pop2 = pio.read(poppath) # this should be the same population as pop

# pass a population to the model to use that instead of creating a new one
model3 = TITAN(param2, pop2)
model3.run(outdir)
```

### Running the Tests

To make sure everything is working, run the tests.  A few of the tests (marked `integration_stochastic`) sometimes fail as they are testing for general behavior and not specific correctness, but all other tests should always pass.

`python -m pytest`
