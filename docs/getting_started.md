# Getting Started

To get started, install the requirements listed in the `requirements.txt` using a local python install or virtual env. Once installed, the model can be run using the `run_titan.py` program and configured using the `params.py` file.

!!! tip
    Running a large job locally? Look into using [pypy](https://www.pypy.org/) instead of python for MOAR performance.  This is what we use on OSCAR.  Otherwise, all of the instructions hold, just using pypy and pypy's pip.

## Installing the Dependencies

In order for TITAN to run, you must install `Python3.x` and the necessary python packages listed in the `requirements.txt`. This can be done by using `pip` via:

```
pip install -r requirements.txt
```

### Dev Dependencies

If you are planning to contribute to the package, also install the dev dependencies.  These packages help with code style, linting, typing, and building documentation.

```
pip install -r dev_requirements.txt
```
