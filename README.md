# TITAN Simulation
[![DOI](https://zenodo.org/badge/80315242.svg)](https://zenodo.org/badge/latestdoi/80315242)
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/pph-collective/TITAN)](https://github.com/pph-collective/TITAN/releases/latest/) [![](https://github.com/pph-collective/TITAN/workflows/Unit%20Tests/badge.svg)](https://github.com/pph-collective/TITAN/actions) [![codecov](https://codecov.io/gh/pph-collective/TITAN/branch/develop/graph/badge.svg?token=wjkExshhyh)](https://codecov.io/gh/pph-collective/TITAN) [![GitHub](https://img.shields.io/github/license/pph-collective/TITAN)](https://github.com/pph-collective/TITAN/blob/develop/LICENSE) [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://pph-collective.github.io/TITAN/)

TITAN (Treatment of Infectious Transmissions through Agent-based Network) is an agent-based simulation model used to explore contact transmission in complex social networks. Starting with the initializing agent population, TITAN iterates over a series of stochastic interactions where agents can interact with one another, transmit infections through various medium, and enter and exit the care continuum. The purpose of TITAN is to evaluate the impact of prevention and treatment models on incidence and prevalence rates of the targeted disease(s) through the use of data fitting simulated trajectories and rich statistics of primary/sub-population attributable proportions.

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

To run the model, execute the `run_titan.py` program within the `/titan/` directory. See [TITAN params](https://pph-collective.github.io/titan-params-app) for documentation on how to set and use parameters.

Results of the model are generated and aggregated into the `/results/` directory by default. If the model is re-run, the existing results will be overwritten. A helper script has been written to prepare simulations for use with OSCAR, and is labelled `subTitan.sh` in the root directory.


### Running the tests

`python -m pytest`

Code coverage is tracked via CodeCov and targets are set for unit, integration, and the overall project.  Any changes to the code should also have associated tests.


## Built With
* [Python3.x](https://www.python.org/downloads/release/python-374/) - Programming language

  Van Rossum G, Drake FL. Python 3 Reference Manual. Scotts Valley, CA: CreateSpace; 2009.

* [Networkx](https://networkx.github.io/) - Network structure backend

  Hagberg A, Swart P, S Chult D. Exploring network structure, dynamics, and function using NetworkX. 2008.

* [Numpy](http://www.numpy.org/) - Numerical libraries

  Oliphant TE. A guide to NumPy. Vol. 1. Trelgol Publishing USA; 2006.

* [Matplotlib]()

  Hunter JD. Matplotlib: A 2D graphics environment. Computing in science &amp; engineering. 2007;9(3):90–5.


## Authors

* **Lars Seeman** - *Initial work*
* **Max King** - *Continued development*
* **Sarah Bessey** - *Continued development*
* **Mary McGrath** - *Continued development*

## License

This project is licensed under the GNU General Public License version 3 - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
