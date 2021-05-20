Locations can be used in TITAN to differentiate agents by "geography".  The primary features of locations are:

* Differentiated parameters via location scaling/overrides
    - Allows different demographics or interventions by location
    - See [params app for details](https://pph-collective.github.io/titan-params-app/#/params#location-1)
* Location based assorting (including based on neighboring locations)
    - Can have agents assort with agents from their own location vs neighbors vs all others
    - Neighboring locations are determined by the `edges` defined in `params.location`
    - See [params app for details on assorting rules](https://pph-collective.github.io/titan-params-app/#/params#assort_mix-1)
    - It is also possible to define edges via a geography CSV (see [utils](/api/utilities/#titan.utils.grid_file_to_edges)). This is also exposed via the `grid2edges` command line utility (run `grid2edges --help` for usage).
* Migration between locations
    - When an agent migrates locations, they adopt the parameters of their new location
    - Migration can cause the population numbers in a location to drift over time
    - See [params app for details](https://pph-collective.github.io/titan-params-app/#/params#location-1)

## Location

::: titan.location.Location

## LocationEdge

::: titan.location.LocationEdge

## Geography

::: titan.location.Geography
