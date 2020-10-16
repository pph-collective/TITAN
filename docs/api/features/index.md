## Agent Features

Features from [`params.features`](https://marshall-lab.github.io/titan-params-app/#/params#features-1) that include an agent attribute (e.g. PrEP, incarceration) are implemented in standalone files that implement the interface from `features.BaseFeature`.  This allows for all of the logic related to a feature to be consolidated in one place and additionally makes incorporating a new feature into the model as simple as possible.

To add a new feature:

* Add the param to `param.features`
* Add a file to `features/` which creates a class that is a sub-class of `BaseFeature`
    * Implement the methods of `BaseFeature` which are needed for this feature
    * Not all methods are needed for all features (see below for details on the methods)
* Re-export the feature from `features/__init__.py`
* Add it to the docs in `docs/api/features/` and to the nav in `mkdocs.yml`

The `HIVModel`, `Population`, and `Agent` classes all use sub-classes of `BaseFeature` to initialize the object/call methods as appropriate.



::: titan.features.base_feature
