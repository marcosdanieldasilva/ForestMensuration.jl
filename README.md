# ForestMensuration

ForestMensuration.jl provides advanced functions for dendrometric calculations in Julia.
Its focus is on accurate **tree cubage (volume estimation)**, **dendrometric averages**
(mean, quadratic mean, dominant diameter/height, ...), and **frequency/diametric
distribution tables**, built on top of
[ForestCore.jl](https://github.com/JuliaForests/ForestCore.jl) so every measurement
carries its unit and is converted automatically. These methods are essential for forest
mensuration and biometrics, supporting forest inventory, management, and research.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://marcosdanieldasilva.github.io/ForestMensuration.jl/stable/forestmensuration/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://marcosdanieldasilva.github.io/ForestMensuration.jl/dev/forestmensuration/)
[![Build Status](https://github.com/marcosdanieldasilva/ForestMensuration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/marcosdanieldasilva/ForestMensuration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/marcosdanieldasilva/ForestMensuration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/marcosdanieldasilva/ForestMensuration.jl)

## Installation

Install the package via Julia's package manager:

```julia-repl
using Pkg
Pkg.add("ForestMensuration")
```

## Overview

ForestMensuration.jl is designed for professionals and researchers in forestry,
dendrometry, and forest biometrics. Its key features include:

- **Cubage (Tree Volume Estimation):**
  Rigorous section methods (Smalian, Huber, Newton, Hohenadl) for individual logs and
  full bole profiles, plus a high-level `cubage` function that partitions a tree into
  commercial, residual, stump, and tip volumes and computes form factors.

- **Dendrometric Averages:**
  Arithmetic and quadratic mean diameter, Weise diameter, central basal-area diameter,
  Hohenadl's diameters, dominant diameter/height (Assmann), and Lorey's mean height —
  individually or bundled into per-plot summary tables (`dmetrics`, `hmetrics`,
  `standmetrics`).

- **Frequency and Diametric Distribution Tables:**
  Class-based frequency tables for any variable (`frequencytable`), and a
  diameter-specific version that also reports basal area and per-hectare/per-acre
  expansions (`diametrictable`).

Every function accepts `Unitful` quantities directly, and also accepts plain numbers for
convenience — diameters and bark thickness default to centimeters, heights and log
lengths to meters, volumes to cubic meters, and plot areas to hectares. Results always
carry units.

## Example Usage

### Cubage Calculation

Compute the volume of a single log using the Smalian method, or partition a whole tree
into a detailed volume/form-factor table.

```julia-repl
using ForestMensuration

# a single log: 3 m long, 30 cm base diameter, 25 cm top diameter
julia> smalian(3.0u"m", 30.0u"cm", 25.0u"cm")
0.17965982987716633 m^3

# a full tree profile: cumulative heights and diameters from the ground up
julia> h = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3]u"m";
julia> d = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9]u"cm";
julia> cubage(h, d)
1×11 DataFrame
 Row │ vt            v0              vc             vr         vn               d          h          hc         aff       nff       qf
     │ Quantity…     Quantity…       Quantity…      Quantity…  Quantity…        Quantity…  Quantity…  Quantity…  Float64   Float64   Float64
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 0.022122 m^3  0.00190852 m^3  0.0194575 m^3    0.0 m^3  0.000756077 m^3     7.0 cm      9.3 m      7.3 m  0.618097  0.505558  0.761071
```

### Dendrometric Averages

Summarize the diameter and height structure of a plot, individually or per-hectare.

```julia-repl
using ForestMensuration

julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]u"cm";
julia> heights = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2]u"m";
julia> plotarea = 0.05u"ha";

julia> dmetrics(diameters, plotarea)
1×8 DataFrame
 Row │ dl          dm         dg          dw         dz          dd         du          dv
     │ Quantity…   Quantity…  Quantity…   Quantity…  Quantity…   Quantity…  Quantity…   Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────
   1 │ 12.7085 cm   17.25 cm  17.7799 cm    18.6 cm  17.2663 cm    21.0 cm  21.7915 cm  26.3274

# a full plot summary: absolute totals, per-hectare density/basal area, and every diameter/height metric
julia> standmetrics(diameters, heights, plotarea)
1×19 DataFrame
 Row │ n      g             EF          N            G                  dl          dm         dg          dw         dz          dd         du          dv       hl         hm         hd         hg         hu         hv
     │ Int64  Quantity…     Quantity…   Quantity…    Quantity…          Quantity…   Quantity…  Quantity…   Quantity…  Quantity…   Quantity…  Quantity…   Float64  Quantity…  Quantity…  Quantity…  Quantity…  Quantity…  Float64
─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │    10  0.248284 m^2  20.0 ha^-1  200.0 ha^-1  4.96568 m^2 ha^-1  12.7085 cm   17.25 cm  17.7799 cm    18.6 cm  17.2663 cm    21.0 cm  21.7915 cm  26.3274  11.9589 m    15.55 m     18.5 m   17.148 m  19.1411 m   23.094
```

### Frequency and Diametric Tables

Build a class-based frequency table for any variable, or a diameter-specific table with
basal area and per-hectare expansion.

```julia-repl
using ForestMensuration

julia> frequencytable(diameters, 2u"cm")
8×7 DataFrame
 Row │ LI         Xi         LS         fi     Fi     fri      Fri
     │ Quantity…  Quantity…  Quantity…  Int64  Int64  Float64  Float64
─────┼─────────────────────────────────────────────────────────────────
   1 │   10.0 cm    11.0 cm    12.0 cm      1      1     10.0     10.0
   2 │   12.0 cm    13.0 cm    14.0 cm      1      2     10.0     20.0
   3 │   14.0 cm    15.0 cm    16.0 cm      1      3     10.0     30.0
   4 │   16.0 cm    17.0 cm    18.0 cm      2      5     20.0     50.0
   5 │   18.0 cm    19.0 cm    20.0 cm      1      6     10.0     60.0
   6 │   20.0 cm    21.0 cm    22.0 cm      2      8     20.0     80.0
   7 │   22.0 cm    23.0 cm    24.0 cm      1      9     10.0     90.0
   8 │   24.0 cm    25.0 cm    26.0 cm      1     10     10.0    100.0

julia> diametrictable(diameters, 2u"cm", plot_area=plotarea)
8×14 DataFrame
 Row │ LI         Xi         LS         fi     Fi     fri      Fri      g               ng              ∑ng             fi_ha       Fi_ha        ng_ha               ∑ng_ha
     │ Quantity…  Quantity…  Quantity…  Int64  Int64  Float64  Float64  Quantity…       Quantity…       Quantity…       Quantity…   Quantity…    Quantity…           Quantity…
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │   10.0 cm    11.0 cm    12.0 cm      1      1     10.0     10.0  0.00950332 m^2  0.00950332 m^2  0.00950332 m^2  20.0 ha^-1   20.0 ha^-1  0.190066 m^2 ha^-1  0.190066 m^2 ha^-1
   2 │   12.0 cm    13.0 cm    14.0 cm      1      2     10.0     20.0   0.0132732 m^2   0.0132732 m^2   0.0227765 m^2  20.0 ha^-1   40.0 ha^-1  0.265465 m^2 ha^-1  0.455531 m^2 ha^-1
  ⋮  │     ⋮          ⋮          ⋮        ⋮      ⋮       ⋮        ⋮            ⋮               ⋮               ⋮             ⋮           ⋮               ⋮                   ⋮
```

## Keywords

Forest mensuration, dendrometry, forest inventory, tree cubage, dendrometric averages,
frequency distribution, diametric distribution, forest biometrics.

## License

This project is licensed under the MIT License.
