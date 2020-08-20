## Running the Model

To run the model, execute the `run_titan.py` program within the `/titan/` directory. See [TITAN params](https://marshall-lab.github.io/titan-params-app) for documentation on how to set and use parameters.

Results of the model are generated and aggregated into the `/results/` directory by default. If the model is re-run, the existing results will be overwritten. A helper script has been written to prepare simulations for use with OSCAR, and is labelled `subTitan.sh` in the root directory.


### Running the tests

`python -m pytest`

Code coverage is tracked via CodeCov and targets are set for unit, integration, and the overall project.  Any changes to the code should also have associated tests.
