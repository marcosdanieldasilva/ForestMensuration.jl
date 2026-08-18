"""
Description:

ForestMensuration.jl is a Julia package that provides a comprehensive set of functions for performing dendrometric calculations. Designed with ease of use in mind, the package offers tools that simplify complex forestry calculations, making it straightforward to:

- **Calculate tree and stand volume (cubage)**: Support for various methods such as Huber, Smalian, and Newton allows precise calculation of tree and stand volumes.
- **Compute dendrometric averages**: Calculate essential dendrometric metrics like mean diameter, quadratic mean diameter, and others to understand stand structure.
- **Create frequency tables**: Generate frequency and diametric tables to analyze the distribution of dendrometric variables such as diameter and height.
- **Estimate forest inventory sampling designs**: All 10 classic designs — simple
  (`simplecasualsampling`), stratified (`stratifiedsampling`), systematic
  (`systematicsampling`, `multistartsystematicsampling`), cluster/two-stage
  (`clustersampling`, `twostagesampling`), and the four sampling-on-successive-occasions
  designs (`independentoccasionssampling`, `completereplacementsampling`,
  `partialreplacementsampling`, `doublesampling`) — every one generic to any number of
  strata/clusters/plots and returning a [`SamplingReport`](@ref) or plain `DataFrame`.
- **Apply fitted stem taper models**: [`taperdiameter`](@ref)/[`taperheight`](@ref) evaluate
  or invert a [`TaperFit`](@ref) (fit with `ForestModeling.fit` — see its module docstring
  for the 10-model catalog), [`taperedvolume`](@ref) integrates volume from the fitted
  curve, and [`logassortment`](@ref) simulates cutting a stem into logs across any number
  of product classes.
- **Classify site productivity**: [`siteClassification`](@ref)/[`hdomClassification`](@ref)/
  [`siteTable`](@ref) — guide-curve (anamorphic, delta-method) site index classification
  from any allometric model fitted with age as its sole continuous regressor.
  `siteTable`'s automatic class breadth reuses the same Sturges'-rule-based "nice number"
  rounding as [`frequencytable`](@ref)/[`diametrictable`](@ref).
- **Fit and select allometric regression models**: Re-exported from
  [ForestModeling.jl](https://github.com/JuliaForests/ForestModeling.jl) — search a bounded
  catalog of hypsometric/volumetric transforms (`regression`), rank the results
  (`criteriaTable`/`criteriaSelection`), and compare pooled/grouped/covariate-adjusted
  variants (`regressionGrouped`), with range-safe prediction across them
  (`predictBounded`).

The package facilitates the analysis of dendrometric and forest data, performs complex calculations with simple commands, and offers a user-friendly and intuitive interface.

# Units of measurement

Every function is built on `Unitful` quantities, so measurements carry their unit and are
converted automatically. Diameters may be given in `cm`, `mm`, `inch`, ... and heights in
`m`, `ft`, ..., in any combination; results are normalised to a single coherent unit:
areas come back in `m^2` and volumes in `m^3`, or in `ft^2`/`ft^3` when the diameters are
imperial.

For convenience every function also accepts plain numbers, in which case the package
assumes its default units:

| Measurement                 | Assumed unit |
|:----------------------------|:-------------|
| diameter, bark thickness    | `cm`         |
| height, log length          | `m`          |
| volume                      | `m^3`        |
| plot area                   | `ha`         |

So `cubage([0.3, 1.3, 3.3], [9.0, 7.0, 5.8])` reads the heights as meters and the
diameters as centimeters, and returns exactly the same table as the equivalent call with
explicit units. The one deliberate exception is `frequencytable`, which classifies any
variable and therefore stays unitless when given unitless data — use `diametrictable`
for diameter-specific tables.

Results always carry units; use `removeunits`/`restoreunits` (re-exported from
`ForestFoundations`) to move between unitful DataFrames and plain numeric ones for export.
"""
module ForestMensuration
using DataFrames, Distributions, ForestFoundations, ForestModeling, LinearAlgebra, QuadGK, Reexport, Roots, StatsBase, Tables

import StatsModels: modelcols, FunctionTerm

@reexport using ForestFoundations
@reexport using ForestModeling

include("distributiontables.jl")
include("siteclassification.jl")
include("dendrometrics.jl")
include("cubage.jl")
include("inventory/common.jl")
include("inventory/inventoryreport.jl")
include("inventory/simplecasualsampling.jl")
include("inventory/stratifiedsampling.jl")
include("inventory/systematicsampling.jl")
include("inventory/clustersampling.jl")
include("inventory/multistartsystematicsampling.jl")
include("inventory/twostagesampling.jl")
include("inventory/independentoccasionssampling.jl")
include("inventory/completereplacementsampling.jl")
include("inventory/partialreplacementsampling.jl")
include("inventory/doublesampling.jl")
include("taper/common.jl")
include("taper/taperdiameter.jl")
include("taper/taperheight.jl")
include("taper/taperedvolume.jl")
include("taper/logassortment.jl")

export
  # Cubage
  artificialformfactor,
  barkfactor,
  conevolume,
  cubage,
  cylindervolume,
  diameterinterpolation,
  heightinterpolation,
  naturalformfactor,
  quotientform,
  smalian,
  newton,
  huber,
  hohenadl,
  # Frequency and Statistic functions
  dmetrics,
  hmetrics,
  standmetrics,
  diametrictable,
  frequencytable,
  # Forest inventory sampling
  SamplingReport,
  simplecasualsampling,
  stratifiedsampling,
  systematicsampling,
  clustersampling,
  multistartsystematicsampling,
  twostagesampling,
  independentoccasionssampling,
  completereplacementsampling,
  partialreplacementsampling,
  doublesampling,
  # Stem taper application (fitting lives in ForestModeling.jl, re-exported above)
  taperdiameter,
  taperheight,
  taperedvolume,
  logassortment,
  # Site classification
  siteClassification,
  hdomClassification,
  siteTable

end