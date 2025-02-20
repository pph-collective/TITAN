We welcome issues and pull requests to help improve the TITAN model.  We use the following guidelines to ensure standards are met.

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

### pypy

pypy is a JIT compiled version of python which generally makes code faster.  It is important that titan remain installable/usable with pypy and we test for this with GitHub actions.  However, pypy doesn't play well with some of our linting tools, so we don't typically use it for development.


## GitHub Workflow

When working on TITAN, make a branch from `develop` to make changes in, then make a pull request to `develop` when ready for review.

### Branches

* `develop` is the primary working branch for the project
* `main` is the branch releases are made off of and is in stable condition at all times
* topic branches are created for new features, fixes, or really any changes

### Commitizen

All commits should use the angular style.  Commitizen makes this easy.

```
pip install --user commitizen
```

Then use `cz commit` where you would have previously done `git commit`.

## Code Standards

### Testing

We strive to have test coverage for all of the features of TITAN.  When adding a new feature or fixing a bug, add tests for the new feature or that test the bug condition to make sure that the bug doesn't crop up agin.

[Pytest](https://docs.pytest.org/en/stable/) is our test runner.  It runs all of our unit and integration tests and reports back any errors or failures.  These should almost always pass (the almost refers to our stochastic integration tests which have some element of randomness in their success).

```
poetry run pytest # runs all tests
poetry run pytest -m unit # only unit tests
poetry run pytest -m integration_deterministic # only deterministic integration tests
poetry run pytest -m integration_stochastic # only stochastic integration tests
```

### Code Style

### black

The code must conform to `black`'s standards and this is automatically checked via github actions.  To automatically fix files, run `black .` from the root of the `TITAN` directory.

### flake8
The code must not have any egregious linting errors. And others should be minimized as reasonable.

Check for egregious errors:
```
flake8 titan --count --select=E9,F63,F7,F82 --show-source --statistics
```

Check for all warnings:
```
flake8 titan --count --exit-zero --max-complexity=12 --max-line-length=88 --statistics
```

### Typing

Please use [type hints](https://mypy.readthedocs.io/en/stable/) on all signatures where reasonable.  This will make sure the code is more readable, can be statically tested for type soundness, and helps fill in the documentation.

Run the below to check for type errors:
```
mypy titan
```

### Documentation

All functions and methods should have an up to date [google style docstring](https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html).  These docstrings are used to build TITAN's documentation website.  Additional prose can be added to the website by editing the appropriate markdown file in the `docs/` directory.

To develop/see the docs locally, run:
```
mkdocs serve
```
then navigate to http://localhost:8000 to see your docs.  They will automoatically reload as you make changes to files.
