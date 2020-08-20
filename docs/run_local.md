## Running the Model

The model has a wrapper script called `run_titan.py` that makes running a full simulation easy.  TITAN can also be run from an interactive repl or a custom script.

### run_titan.py

To run the model, execute the `run_titan.py` program within the `/titan/` directory. See [TITAN params](https://marshall-lab.github.io/titan-params-app) for documentation on how to set and use parameters.

Results of the model are generated and aggregated into the `/results/` directory by default. If the model is re-run, the existing results will be overwritten.

#### Usage

Below are the results of `python run_titan.py --help`.  It highlights all of the command line arguments that can be passed to the script.

```
usage: run_titan.py [-h] [-n [NMC]] [-S SETTING] -p PARAMS [-o OUTDIR]
                    [-b BASE] [-e] [--savepop SAVEPOP] [--poppath POPPATH]
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
  --savepop SAVEPOP     Save population after creation, but before model run.
                        'all' = save all atributes, 'core' = save core (non-
                        intervention) attributes.
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

::: run_titan.main
    rendering:
        show_root_heading: true
        show_root_toc_entry: false
        heading_level: 4


### Running Interactively



### Running the Tests

To make sure everything is working, run the tests.

`python -m pytest`
