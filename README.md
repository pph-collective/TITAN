# TITAN Simulation

![](https://github.com/marshall-lab/TITAN/workflows/Unit%20Tests/badge.svg) [![codecov](https://codecov.io/gh/marshall-lab/TITAN/branch/develop/graph/badge.svg?token=wjkExshhyh)](https://codecov.io/gh/marshall-lab/TITAN)

TITAN (Treatment of Infectious Transmissions through Agent-based Network) is an agent-based simulation model used to explore contact transmission in complex social networks. Starting with the initializing agent population, TITAN iterates over a series of stochastic interactions where agents can engage in risk behavior with one another, transmit infections through various medium, and enter and exit the care continuum. The purpose of TITAN is to evaluate the impact of treatment models on incidence and prevalence rates of the targeted disease(s) through the use of data fitting simulated trajectories and rich statistics of primary/sub-population attributable proportions.

Agent populations are defined as graphs (nodes connected by edges). Nodes in the graph are used to represent the attributes (or collection of attributes) of an agent (person), and edges define the type of relationship between agents. In practice, a graph represents a social network of connected people through various relationship types, and provides the medium for which agents can interact.

## Getting Started

To get started, install the requirements listed in the `requirements.txt` using a local python install or virtual env. Once installed, the model can be run using the `run_titan.py` program and configured using the `params.py` file.

### Prerequisites

In order for TITAN to run, you must install `Python3.x` and the necessary python packages listed in the `requirements.txt`. This can be done by using `pip` via:

```
pip install -r requirements.txt
```

### Installing

Currently, TITAN does not require any further installation as the source code exists within this repo in `/titan`

## Running the Model

To run the model, execute the `run_titan.py` program within the `/titan/` directory. It will natively read in the `params.py` configuration file as the input parameters for the model. These can be configured as desired for the target model setting and configuration. Other examples of calibrated settings exist within the `/settings/` directory.

Results of the model are generated and aggregated into the `/results/` directory. If the model is re-run, the existing results will be overwriten. A helper script has been written to prepare simulations for use with OSCAR, and is labelled `subTitan.sh` in the root directory.


### Running the tests

`python -m pytest`

Code coverage is tracked via CodeCov and targets are set for unit, integration, and the overall project.  Any changes to the code should also have associated tests.


## Built With
* [Python3.x](https://www.python.org/downloads/release/python-374/) - Programming language
* [Networkx](https://networkx.github.io/) - Network structure backend
* [Numpy](http://www.numpy.org/) - Numerical libraries


## Authors

* **Lars Seeman** - *Initial work*
* **Max King** - *Continued development*
* **Sarah Bessey** - *Continued development*
* **Mary McGrath** - *Continued development*

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
