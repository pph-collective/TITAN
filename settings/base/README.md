# Base Settings

These parameters are more complex defaults than can be expressed via a param definition file.  Primarily these are cases where there are established defaults that differ across subsets of the population.

`create_params` and `run_titan` use these settings by default, but they can be turned off by passing `use_base=False`.

#### Order of Param Precedence

1. Param
2. Setting
3. Base
4. Param Definition
