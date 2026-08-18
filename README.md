# ForestMensuration

ForestMensuration.jl provides advanced functions for dendrometric calculations in Julia.
Its focus is on accurate **tree cubage (volume estimation)**, **dendrometric averages**
(mean, quadratic mean, dominant diameter/height, ...), **frequency/diametric
distribution tables**, **forest inventory sampling** (all 11 classic designs), and
**allometric regression/site classification**, built on top of
[ForestFoundations.jl](https://github.com/JuliaForests/ForestFoundations.jl) (units) and
[ForestModeling.jl](https://github.com/JuliaForests/ForestModeling.jl) (the regression
engine, re-exported here) so every measurement carries its unit and is converted
automatically. These methods are essential for forest mensuration and biometrics,
supporting forest inventory, management, and research.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaForests.github.io/ForestMensuration.jl/stable/forestmensuration/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaForests.github.io/ForestMensuration.jl/dev/forestmensuration/)
[![Build Status](https://github.com/JuliaForests/ForestMensuration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaForests/ForestMensuration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaForests/ForestMensuration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaForests/ForestMensuration.jl)

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

- **Forest Inventory Sampling:**
  All 11 classic sampling designs, each generic to any number of strata/clusters/plots:
  simple (`simplecasualsampling`), stratified (`stratifiedsampling`), systematic
  (`systematicsampling`, `multistartsystematicsampling`), cluster/two-stage
  (`clustersampling`, `twostagesampling`), horizontal point/Bitterlich angle-count
  sampling (`horizontalpointsampling`), and the four sampling-on-successive-occasions
  designs used in continuous forest inventory (`independentoccasionssampling`,
  `completereplacementsampling`, `partialreplacementsampling`, `doublesampling`).

- **Stem Taper Equations:**
  Fit 10 classic published taper (stem-profile) forms — `Kozak1969`, `Schoepfer1966`,
  `Matte1949`, `Demaerschalk1972`, `Clutter1980`, `MaxBurkhart1976`, `Johnson1911`,
  `Kozak1988`, `Kozak2004`, `Bi2000` — with `fit` (re-exported from
  [ForestModeling.jl](https://github.com/JuliaForests/ForestModeling.jl)), then evaluate
  diameter/height along the stem (`taperdiameter`, `taperheight`), integrate volume between
  two heights (`taperedvolume`), and simulate cutting the stem into logs across any number
  of product classes (`logassortment`).

- **Site Classification:**
  Guide-curve (anamorphic, delta-method) site index classification from any allometric
  model fitted with age as its sole continuous regressor (`siteClassification`,
  `hdomClassification`, `siteTable`) — `siteTable`'s automatic class breadth reuses the
  same rounding convention as `frequencytable`/`diametrictable`.

- **Allometric Regression** (re-exported from
  [ForestModeling.jl](https://github.com/JuliaForests/ForestModeling.jl)):
  Search a bounded catalog of hypsometric/volumetric transforms and rank the results
  (`regression`, `criteriaTable`, `criteriaSelection`), and compare pooled/grouped/
  covariate-adjusted variants (`regressionGrouped`) with range-safe prediction across
  them (`predictBounded`).

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

### Regression / Allometric Modeling

Search a bounded catalog of hypsometric transforms, rank the results, and predict —
units flow through automatically, just like every other function in this package.

```julia-repl
using ForestMensuration

julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]u"cm";
julia> heights = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2]u"m";
julia> data = DataFrame(d=diameters, h=heights);

julia> models = regression(data, :h, :d; nMax=2);

julia> best = criteriaSelection(models, :adjr2, :cv)
h ^ -1 = 0.0582 - 3.338e-5 * d ^ 2 + 4.861 * d ^ -2

julia> predict(best)
10-element Vector{Union{Missing, Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}}:
 10.147829089308829 m
 11.483629283168705 m
                   ⋮
 21.130327937569657 m
```

#### Robust Regression

For heavy outliers or strongly relative error, `fitRobust` minimizes an alternative loss
(`SSE`, `MAE`, `HUBER`, `MSLE`, `MAPE`) with `Optim.jl` instead of closed-form least
squares — units flow through the same way:

```julia-repl
julia> mr = fitRobust(@formula(log(h) ~ log(d)), data, HUBER)
log(h) = 0.3021 + 0.8593 * log(d)

julia> predict(mr)
10-element Vector{Union{Missing, Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}}:
 10.204432130945504 m
 11.445132745669772 m
                   ⋮
 20.76315659021014 m
```

#### Grouped (Stratified) Regression

For stratified data (e.g. several species that plausibly need different equations),
`regressionGrouped` compares one pooled equation, one pooled equation with the group as a
covariate, and fully separate equations per group:

```julia-repl
julia> gdata = DataFrame(d=..., h=..., species=...);  # 3 species, distinct slope/intercept each

julia> gm = regressionGrouped(gdata, :h, :d, :species);

julia> criteriaTable(gm, :adjr2, :cv)
3×4 DataFrame
 Row │ type     model                              adjr2     cv
     │ String   Allometr…?                         Float64   Float64
─────┼────────────────────────────────────────────────────────────────
   1 │ General  h = 0.3539 - 0.03555 * d + 0.660…  0.822337  5.38691
   2 │ Qualy    h ^ -1 = 0.8596 + 0.05805 * √d -…  0.995854  0.822871
   3 │ Grouped  missing                            0.997638  0.621154
```

Here stratifying by species clearly pays off: `Grouped` and `Qualy` both dominate
`General`. `predictBounded` scores new data across the three per row, falling back to the
broader model whenever a row's predictor falls outside the range its own group's equation
was actually fit on — so a small subgroup can't extrapolate into nonsensical predictions:

```julia-repl
julia> predictBounded(gm, DataFrame(d=[15.0, 60.0, 15.0], species=["Cedar", "Cedar", "Birch"]))
```

See the
[ForestModeling.jl README](https://github.com/JuliaForests/ForestModeling.jl#readme) for
the full regression API, including named classic equations (`fitClassic`) and inference
accessors (`coef`, `confint`, `cooksdistance`, ...).

### Site Classification

Guide-curve (anamorphic, delta-method) site index classification from any allometric
model fitted with age as its sole continuous regressor — `siteClassification` gives the
site index at a chosen index age, `hdomClassification` inverts it (site index → expected
dominant height), and `siteTable` tabulates dominant height by age × site class.
`siteTable`'s automatic class breadth reuses the same "nice number" Sturges'-rule rounding
as [`frequencytable`](#frequency-and-diametric-tables).

```julia-repl
using ForestMensuration, DataFrames

julia> data = DataFrame(
         plot = repeat(1:6, inner=5),
         age  = repeat([36, 48, 60, 72, 84], outer=6),
         h    = [13.6, 17.8, 21.5, 21.5, 21.8, 14.3, 17.8, 21.0, 21.0, 21.4,
                 14.0, 17.5, 21.2, 21.2, 21.4, 13.4, 18.0, 20.8, 20.8, 23.2,
                 13.2, 17.4, 20.3, 20.3, 22.0, 13.2, 17.8, 21.3, 21.3, 22.5]u"m",
       );

julia> reg = criteriaSelection(regression(data, :h, :age), :adjr2, :cv)
h ^ -1 = -0.0424 + 0.01769 * log(age) + 68.13 * age ^ -2

julia> site = siteClassification(reg, data, 60)   # site index at index age 60, carries h's unit (m)
30-element Vector{Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}:
 20.5 m
 20.2 m
  ⋮
 21.0 m

julia> siteTable(reg, 60)   # predicted dominant height per age × site class; stays plain Float64 (meters) even here
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

### Forest Inventory Sampling

Estimate a stand's total volume — with its confidence interval, coefficient of variation,
and required sample size — from any of the 11 classic sampling designs. Every design is
generic to any number of strata/clusters/plots and accepts plain numbers (volume defaults
to `m^3`, areas to `ha`) or explicit `Unitful` quantities, exactly like the rest of the
package.

Every result table is a single-row, one-column-per-statistic `DataFrame` — the same shape
as `dmetrics`/`standmetrics` — so every value that carries a physical unit is a real
`Unitful` quantity rather than a plain number paired with a separate units column. That
means `removeunits`/`restoreunits` (re-exported from ForestFoundations.jl) work on it
directly, just like on any other table in this package.

```julia-repl
using ForestMensuration

julia> v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1];

julia> report = simplecasualsampling(v, 0.05, 10; e=10, α=0.95)
1×17 DataFrame
 Row │ vm           cv       s2m          se           abserr       relerr   vha                vtotal       cilower      ciupper      pop      f        n      nreq   nmiss  N      area
     │ Quantity…    Float64  Quantity…    Quantity…    Quantity…    Float64  Quantity…          Quantity…    Quantity…    Quantity…    String   Float64  Int64  Int64  Int64  Int64  Quantity…
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 441.691 m^3    10.03  168.498 m^6  12.9807 m^3  28.9227 m^3     6.55  8833.82 m^3 ha^-1   88338.2 m^3  82553.6 m^3  94122.7 m^3  finite    0.945     11      7      0    200  10 ha

julia> report.vm
441.6909090909091 m^3

julia> removeunits(report)[:, ["vm (m^3)", "vtotal (m^3)"]]
1×2 DataFrame
 Row │ vm (m^3)  vtotal (m^3)
     │ Float64   Float64
─────┼────────────────────────
   1 │  441.691       88338.2
```

Multi-table designs (stratified, cluster, two-stage, and the four successive-occasions
designs) return a `SamplingReport` instead of a bare `DataFrame`, so the per-group/
per-occasion tables travel together — each one reachable as a property:

```julia-repl
using DataFrames

julia> data = DataFrame(
         stratum=[1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3],
         volume=[18.2, 21.4, 19.8, 20.1, 32.5, 35.1, 30.8, 12.4, 11.9, 13.6, 12.8, 13.1],
       );

julia> report = stratifiedsampling(:stratum, :volume, 0.1, [12.0, 8.0, 20.0], data);

julia> resultTable(report).vm
18.9025 m^3

julia> resultTable(report).nh   # per-stratum measured plots, as a tuple
(4, 3, 5)

julia> auxiliaryTable(report)   # per-stratum n, mean, variance, and allocation weights
julia> anova(report)             # test for a difference between strata means
```

The remaining nine designs follow the same pattern — see the
[Forest Inventory Sampling tutorial](https://JuliaForests.github.io/ForestMensuration.jl/dev/tutorial/#Forest-Inventory-Sampling)
for a complete, runnable example of every one:

| Function | Design |
|:---------|:-------|
| `simplecasualsampling` | Simple random sampling |
| `stratifiedsampling` | Stratified random sampling |
| `systematicsampling` | Systematic sampling (method of successive differences) |
| `multistartsystematicsampling` | Systematic sampling with multiple random starts |
| `clustersampling` | One-stage cluster sampling |
| `twostagesampling` | Two-stage sampling |
| `horizontalpointsampling` | Horizontal point sampling (Bitterlich's angle-count method) |
| `independentoccasionssampling` | Successive occasions, independent samples |
| `completereplacementsampling` | Successive occasions, complete replacement (matched plots) |
| `partialreplacementsampling` | Successive occasions, partial replacement |
| `doublesampling` | Successive occasions, double sampling with regression |

### Stem Taper Equations

Fit a stem taper curve to measured `(height, diameter)` pairs — one or more trees, each
row is one measured section — then evaluate diameter/height along the stem, integrate
volume, or simulate cutting the stem into logs. `fit` (re-exported from
[ForestModeling.jl](https://github.com/JuliaForests/ForestModeling.jl)) accepts any of 10
classic published taper forms; everything downstream is native to this package.

```julia-repl
using ForestMensuration

julia> dbh = [20.0, 20.0, 20.0, 30.0, 30.0, 30.0, 25.0, 25.0, 25.0, 35.0, 35.0, 35.0];

julia> height = [18.0, 18.0, 18.0, 22.0, 22.0, 22.0, 20.0, 20.0, 20.0, 24.0, 24.0, 24.0];

julia> hi = [0.3, 6.0, 14.0, 0.3, 8.0, 18.0, 0.3, 7.0, 16.0, 0.3, 9.0, 20.0];

julia> di = [22.4, 15.8, 6.1, 33.6, 24.2, 8.7, 27.9, 18.6, 7.4, 38.9, 27.1, 9.9];

julia> ft = fit(Kozak1969(), dbh, height, hi, di);

julia> r2(ft), adjr2(ft)
(0.9968593430987872, 0.9961614193429621)

julia> taperdiameter(ft, 25.0u"cm", 20.0u"m", 7.0u"m")   # diameter at 7 m up a 25 cm dbh, 20 m tree
19.575970870808078 cm

julia> taperheight(ft, 25.0u"cm", 20.0u"m", 15.0u"cm")   # inverse: height at a 15 cm diameter
10.514018504633748 m

julia> taperedvolume(ft, 25.0, 20.0)   # whole-stem volume, integrated from the fitted curve
0.4675569019667837 m^3
```

[`logassortment`](https://JuliaForests.github.io/ForestMensuration.jl/dev/reference/#ForestMensuration.logassortment)
simulates cutting the stem into logs across any number of product classes, cut in table
order (most valuable first) until each product's minimum small-end diameter is reached — a
generalized, bug-fixed equivalent of the R package `timbeR`'s per-model bucking functions:

```julia-repl
julia> using DataFrames

julia> products = DataFrame(
         name=["Sawlog", "Pulpwood"], sed=[18.0, 8.0], minlength=[2.5, 2.0], maxlength=[4.0, 3.0], kerf=[0.03, 0.03],
       );

julia> logassortment(ft, 30.0u"cm", 22.0u"m", products)
1×5 DataFrame
 Row │ product                 volume                 logs    totalvolume   totallogs
     │ Tuple…                  Tuple…                 Tuple…  Quantity…     Int64
─────┼────────────────────────────────────────────────────────────────────────────────
   1 │ ("Sawlog", "Pulpwood")  (0.634265, 0.0894935)  (3, 2)  0.723758 m^3          5
```

`volume`/`logs` are `Tuple`-typed columns (one entry per product, in table order) — see the
[tutorial](https://JuliaForests.github.io/ForestMensuration.jl/dev/tutorial/#Stem-Taper-Equations)
for the full walkthrough, including the model catalog and `criteriaTable`-based selection.

## Keywords

Forest mensuration, dendrometry, forest inventory, tree cubage, dendrometric averages,
frequency distribution, diametric distribution, forest inventory sampling, stratified
sampling, cluster sampling, systematic sampling, successive occasions, continuous forest
inventory, allometric regression, grouped regression, site classification, stem taper
equations, log assortment, sortimentos, forest biometrics.

## License

This project is licensed under the MIT License.
