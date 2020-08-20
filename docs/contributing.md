We welcome issues and pull requests to help improve the TITAN model.  We use the following guidelines to ensure standards are met.

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

`pytest` is the library used for testing.

To run all of the tests:
```
python -m pytest
```

To run only the unit tests:
```
python -m pytest -m unit
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
