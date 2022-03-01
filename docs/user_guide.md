# User Guide

## Installation

To get started, install titan using a local python (version 3.6 or later) install or virtual env using `pip`.

```
pip install titan-model
```

!!! tip
    Running a large job locally? Look into using [pypy](https://www.pypy.org/) instead of python for MOAR performance. Otherwise, all of the instructions hold, just using pypy and pypy's pip.

## Development Setup

Planning to help work on the titan codebase? Here are the tools and processes you need to know about:

### Git/Github

We use [git](https://git-scm.com/) for version control and GitHub for our remote repository. To get started, clone the repository to your local machine.

```
git clone https://github.com/pph-collective/TITAN.git
```

We use angular commits to standardize our commits and encourage better messages. [Commitizen](https://pypi.org/project/commitizen/) makes this easy.  Once installed (would recommend doing this globally as opposed to just for this project), run `cz commit` instead of `git commit`.

### Poetry

[Poetry](https://python-poetry.org/) is a python packaging and dependency management tool.  We use this to install/add/remove dependencies, build our package, and publish it to pypy.  Install poetry per [the install instructions](https://python-poetry.org/docs/#installation), then complete the below steps.

```
poetry install -E linting -E docs
```

The `-E` flags tell poetry that we want the optional dependencies that we use for linting (mypy, black, flake8) and for building documentation.

We recommend using a recent python version for development (check out [pyenv](https://github.com/pyenv/pyenv) for a way to manage different python versions).  You may need to tell poetry which installation of python to use - see [their instructions](https://python-poetry.org/docs/managing-environments/).

Poetry installs all of the dependencies into a virtual environment so that they are isolated from other projects you are working on.  To run any shell commands in the environment, prefix them with `poetry run` (e.g. `poetry run run_titan -p my_params.yml` or `poetry run pytest`) or start a poetry shell with `poetry shell`.

### black

[Black](https://github.com/psf/black) formats our files nicely. Our files must be formatted nicely to be merged back into the titan repo.

```
poetry run black .
```

### flake8

[Flake8](https://flake8.pycqa.org/en/latest/) checks for common linting mistake (e.g. unused imports or variables).  We have some checks that are required to be fixed and some that are optional.

```
# stop the build if there are Python syntax errors or undefined names
poetry run flake8 titan --count --select=E9,F63,F7,F82 --show-source --statistics
# exit-zero treats all errors as warnings. The GitHub editor is 127 chars wide
poetry run flake8 titan --count --exit-zero --max-complexity=12 --max-line-length=88 --statistics
```

### mypy

[Mypy](https://mypy.readthedocs.io/en/stable/) is a static type checker.  It will check the [typing](https://docs.python.org/3/library/typing.html) annotations and our code to make sure everything makes sense (e.g. we don't pass a string to a function that requires a number).

```
poetry run mypy titan
```

### pytest

[Pytest](https://docs.pytest.org/en/stable/) is our test runner.  It runs all of our unit and integration tests and reports back any errors or failures.  These should almost always pass (the almost refers to our stochastic integration tests which have some element of randomness in their success).

```
poetry run pytest # runs all tests
poetry run pytest -m unit # only unit tests
poetry run pytest -m integration_deterministic # only deterministic integration tests
poetry run pytest -m integration_stochastic # only stochastic integration tests
```

### pypy

pypy is a JIT compiled version of python which generally makes code faster.  It is important that titan remain installable/usable with pypy and we test for this with GitHub actions.  However, pypy doesn't play well with some of our linting tools, so we don't typically use it for development.
