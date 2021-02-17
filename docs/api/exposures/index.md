## Agent Exposures

Exposures (transmissible items such as HIV and knowledge) from [`params.exposures`](https://pph-collective.github.io/titan-params-app/#/params#exposures-1) are implemented in standalone files that implement the interface from `exposures.BaseExposure`.  This allows for all of the logic related to an exposure to be consolidated in one place and additionally makes incorporating a new exposure into the model as simple as possible.

To add a new exposure:

* Add the param to `param.exposures`
* Add a file to `exposures/` which creates a class that is a sub-class of `BaseExposure`
    * Implement the methods of `BaseExposure` which are needed for this exposure
    * Not all methods are needed for all exposures (see below for details on the methods)
* Re-export the exposure from `exposures/__init__.py`
* Add tests in `tests/exposures/`
* Add it to the docs in `docs/api/exposures/` and to the nav in `mkdocs.yml`

The `TITAN`, `Population`, and `Agent` classes all use sub-classes of `BaseFeature` to initialize the object/call methods as appropriate.


::: titan.exposures.base_exposure
