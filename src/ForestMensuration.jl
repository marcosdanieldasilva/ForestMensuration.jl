"""
Description:

ForestMensuration.jl is a Julia package that provides a comprehensive set of functions for performing dendrometric calculations. Designed with ease of use in mind, the package offers tools that simplify complex forestry calculations, making it straightforward to:

- **Calculate tree and stand volume (cubage)**: Support for various methods such as Huber, Smalian, and Newton allows precise calculation of tree and stand volumes.
- **Compute dendrometric averages**: Calculate essential dendrometric metrics like mean diameter, quadratic mean diameter, and others to understand stand structure.
- **Create frequency tables**: Generate frequency and diametric tables to analyze the distribution of dendrometric variables such as diameter and height.

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
using DataFrames, ForestFoundations, Reexport, StatsBase, Tables, Unitful

import Unitful: Units, Quantity, NoUnits

@reexport using ForestFoundations

include("units.jl")
include("distributiontables.jl")
include("dendrometrics.jl")
include("cubage.jl")

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
  frequencytable

end