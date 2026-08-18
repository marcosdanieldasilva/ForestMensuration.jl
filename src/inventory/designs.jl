"""
    SamplingDesign

Abstract supertype for the forest-inventory sampling designs dispatched through
[`sampling`](@ref) — one empty singleton subtype per design (`SimpleCasualSampling`,
`ClusterSampling`, `HorizontalPointSampling`, ...) below. `subtypes(SamplingDesign)` lists
every design currently available.
"""
abstract type SamplingDesign end

"""
    sampling(design::Type{<:SamplingDesign}, args...; kwargs...)

Single entry point for every forest-inventory sampling design in this package: pass one of
the `SamplingDesign` subtypes below as the first argument to select the design, then that
design's own arguments. Designs are not interchangeable — each keeps its own argument
signature (a plot-level `volume` vector is not the same shape as grouped `cluster`/`volume`
column names, for instance) — so see the chosen type's own docstring (e.g. `?ClusterSampling`)
for its full `Arguments`/`Returns`/`Mathematical basis`/`Examples`.

`methods(sampling)` lists every design × `Real`/`Unitful` overload at once;
`subtypes(SamplingDesign)` lists just the 11 design types:

- [`SimpleCasualSampling`](@ref) — simple random sampling.
- [`StratifiedSampling`](@ref) — stratified random sampling.
- [`SystematicSampling`](@ref) — systematic sampling (method of successive differences).
- [`MultistartSystematicSampling`](@ref) — systematic sampling with multiple random starts.
- [`ClusterSampling`](@ref) — one-stage cluster sampling.
- [`HorizontalPointSampling`](@ref) — horizontal point sampling (Bitterlich's angle-count method).
- [`TwoStageSampling`](@ref) — two-stage sampling.
- [`IndependentOccasionsSampling`](@ref) — successive occasions, independent samples.
- [`CompleteReplacementSampling`](@ref) — successive occasions, complete replacement (matched plots).
- [`PartialReplacementSampling`](@ref) — successive occasions, partial replacement.
- [`DoubleSampling`](@ref) — successive occasions, double sampling with regression.

# Examples
```julia-repl
julia> v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1];

julia> report = sampling(SimpleCasualSampling, v, 0.05, 10; e=10, α=0.95)
julia> resultTable(report).vm
441.6909090909091 m^3
```
"""
function sampling end

"""
    SimpleCasualSampling <: SamplingDesign

    sampling(::Type{SimpleCasualSampling}, volume::AbstractVector{<:Vol}, plot_area::Area, total_area::Area;
             e::Real=10, α::Real=0.95)
    sampling(::Type{SimpleCasualSampling}, volume::AbstractVector{<:Real}, plot_area::Real, total_area::Real; kwargs...)

Performs simple random sampling for forest inventory analysis with specified plot area and total area.

# Description

Calculates various statistical parameters for a forest inventory using simple random
sampling methodology. It is designed to estimate the total volume of timber within a
forest area by analyzing volume measurements from sample plots.

Simple random sampling is a statistical method where each unit (in this case, a plot) has an equal chance of being selected. This method is commonly used in forest inventories to estimate parameters like mean volume per hectare, total volume, and confidence intervals.

# Arguments

- `volume`: vector of volume measurements from each sampled plot, as a `Vol` (`Unitful.Volume`) quantity. Plain numbers are taken to be `u"m^3"`; an explicit unit (e.g. `u"ft^3"`) is preserved through every reported quantity instead.
- `plot_area`: the area of each sample plot, as an `Area` quantity. All plots are assumed to have the same area. A plain number is taken to be hectares.
- `total_area`: the total area of the forest or stand being inventoried, as an `Area` quantity. A plain number is taken to be hectares.
- `e::Real=10`: the desired relative error margin as a percentage (default is 10%). This represents the maximum acceptable error in the estimate relative to the true mean.
- `α::Real=0.95`: the confidence level for the statistical estimates (default is 95%). This determines the confidence interval for the total volume estimate.

# Returns

- `DataFrame`: a single-row report, one column per statistic. Every column that carries a physical unit is a real `Unitful` quantity (exactly like `dmetrics`/`standmetrics` elsewhere in this package), so `removeunits`/`restoreunits` work on it directly. Columns:
  - `vm`: mean plot volume.
  - `cv`: coefficient of variation (%).
  - `s2m`: variance of the mean.
  - `se`: standard error of the mean.
  - `abserr`: absolute sampling error.
  - `relerr`: relative sampling error (%).
  - `vha`: estimated volume per hectare.
  - `vtotal`: estimated total volume.
  - `cilower`/`ciupper`: confidence interval for the total volume.
  - `pop`: `"finite"` or `"infinite"`, the population classification used.
  - `f`: finite-population correction factor.
  - `n`: number of measured plots.
  - `nreq`: required number of plots for the target error `e`.
  - `nmiss`: additional plots still needed (`0` if `n ≥ nreq`).
  - `N`: number of possible plots in the population.
  - `area`: total area.

# Mathematical basis

The mean plot volume and its coefficient of variation:
```math
\\bar{x} = \\frac{\\sum_{i=1}^{n} x_i}{n} \\qquad CV = \\frac{s}{\\bar{x}} \\times 100
```

For an infinite population, the variance of the mean is `s²/n`; for a finite population of
size `N = \\text{total\\_area}/\\text{plot\\_area}`, a finite-population correction `f = 1 - n/N` is applied:
```math
s^2_{\\bar{x}} = \\frac{s^2}{n} \\qquad \\text{or} \\qquad s^2_{\\bar{x}} = \\frac{s^2}{n}\\left(1 - \\frac{n}{N}\\right)
```

The required sample size for the desired error `e` is solved iteratively (see `_requiredsamplesize`), since the t quantile depends on the sample size through its own degrees of freedom.

# Technical description

Simple random sampling is the reference design against which every other sampling method
in this module is compared: stratified, cluster, and systematic sampling all reduce to it
in the degenerate case of a single stratum/cluster/full census. The population is
classified as infinite whenever the sampling fraction is negligible (`f ≥ 0.98`) or the
plot occupies a full reference area (1 ha or 1 ac), in which case the finite-population
correction is dropped from the variance of the mean.

# Examples

```julia-repl
julia> v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1];

julia> report = sampling(SimpleCasualSampling, v, 0.05, 10; e=10, α=0.95)
1×17 DataFrame
 Row │ vm          cv       s2m         se         abserr     relerr   vha             vtotal      cilower     ciupper     pop      f        n      nreq   nmiss  N      area
     │ Quantity…   Float64  Quantity…   Quantity…  Quantity…  Float64  Quantity…       Quantity…   Quantity…   Quantity…   String   Float64  Int64  Int64  Int64  Int64  Quantity…
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 441.691 m^3    10.03  168.498 m^6  12.9807 m^3  28.9227 m^3     6.55  8833.82 m^3 ha^-1  88338.2 m^3  82553.6 m^3  94122.7 m^3  finite    0.945     11      7      0    200  10 ha

julia> report.vm
441.6909090909091 m^3
```
"""
struct SimpleCasualSampling <: SamplingDesign end

"""
    StratifiedSampling <: SamplingDesign

    sampling(::Type{StratifiedSampling}, stratum::Symbol, volume::Symbol, plot_area::Area, strata_area::AbstractVector{<:Area},
             data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    sampling(::Type{StratifiedSampling}, stratum::Symbol, volume::Symbol, plot_area::Real, strata_area::AbstractVector{<:Real},
             data::AbstractDataFrame; kwargs...)

Performs stratified random sampling for forest inventory analysis, generic to any number of strata.

# Description

Stratified random sampling divides the population into `N` non-overlapping subdivisions
(strata) that are each more internally homogeneous than the population as a whole, then
samples independently within each. This usually shrinks the variance of the estimated
mean relative to simple random sampling of the same total size, at the cost of needing
each stratum's area up front.

# Arguments

- `stratum::Symbol`: name of the column identifying each plot's stratum.
- `volume::Symbol`: name of the volume column, as a `Vol` (`Unitful.Volume`) quantity. Plain numbers are taken to be `u"m^3"`.
- `plot_area`: the area of each sample plot, as an `Area` quantity. A plain number is taken to be hectares.
- `strata_area::AbstractVector{<:Area}`: the total area of each stratum, **in the same order as the sorted unique values of the `stratum` column** (ascending). A plain-number vector is taken to be hectares.
- `data::AbstractDataFrame`: the plot-level data.
- `e::Real=10`: desired relative error margin as a percentage (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with three tables: `anova` (test for a difference between
  strata means), `auxiliaryTable` (per-stratum descriptive statistics and allocation
  weights), and `resultTable` (the final stratified estimates, one row, one column per
  statistic — same shape as [`SimpleCasualSampling`](@ref)'s return value). `resultTable`
  columns:
  - `vm`, `cv`, `s2m`, `se`, `abserr`, `relerr`, `vtotal`, `cilower`, `ciupper`, `pop`, `f`: as in [`SimpleCasualSampling`](@ref), using the stratified mean/variance.
  - `nh`, `nreqh`, `nmissh`: per-stratum measured/required/missing plot counts, as tuples in stratum order.
  - `n`, `nreq`: total measured and required plot counts across all strata.
  - `N`: total number of possible plots across all strata.
  - `areah`: each stratum's area, as a tuple of plain numbers in the same unit `strata_area`
    was supplied in (hectares by default). Unlike every other column here, this one is
    **not** `Unitful`-tagged and is left untouched by `removeunits`/`restoreunits` — a
    `Tuple`-typed column has no single `eltype` those functions can strip/restore a unit
    from, so the value is stored already-unitless instead of silently carrying a hidden
    `Unitful` type through the round trip.

# Mathematical basis

The stratified mean and its variance, weighting each stratum by its area share `pₕ = Aₕ/ΣAₕ`:
```math
\\bar{x}_{st} = \\sum_h p_h \\bar{x}_h \\qquad
s^2_{\\bar{x}_{st}} = \\sum_h \\frac{(p_h s_h)^2}{n_h} - \\frac{\\sum_h p_h s_h^2}{N}
```

Optimal (Neyman) allocation distributes the required total sample size across strata in
proportion to `pₕsₕ`, solved iteratively the same way as [`SimpleCasualSampling`](@ref)
(see `_requiredsamplesize`). The confidence interval uses a Satterthwaite-approximated
effective degrees of freedom that accounts for the strata's differing sizes and variances
(Cochran, 1977, eq. 5.35):
```math
df = \\frac{\\left(\\sum_h g_h s_h^2\\right)^2}{\\sum_h \\dfrac{g_h^2 s_h^4}{n_h - 1}}
\\qquad \\text{where} \\qquad g_h = \\frac{N_h(N_h - n_h)}{n_h}
```

# Technical description

With a single stratum this reduces exactly to [`SimpleCasualSampling`](@ref) — the
`anova` table's between-strata line becomes meaningless (0 degrees of freedom) and the
weighted mean collapses to the plain sample mean. The `auxiliaryTable`'s `ps`/`ps²`
columns are the building blocks of every downstream formula (allocation, variance,
required sample size) and are exposed directly so they can be audited.

# Examples

```julia-repl
julia> using DataFrames

julia> data = DataFrame(
         stratum=[1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3],
         volume=[18.2, 21.4, 19.8, 20.1, 32.5, 35.1, 30.8, 12.4, 11.9, 13.6, 12.8, 13.1],
       );

julia> report = sampling(StratifiedSampling, :stratum, :volume, 0.1, [12.0, 8.0, 20.0], data);

julia> resultTable(report)
julia> auxiliaryTable(report)
julia> anova(report)
```
"""
struct StratifiedSampling <: SamplingDesign end

"""
    SystematicSampling <: SamplingDesign

    sampling(::Type{SystematicSampling}, volume::AbstractVector{<:Vol}, plot_area::Area, total_area::Area;
             line::Union{AbstractVector,Nothing}=nothing, e::Real=10, α::Real=0.95)
    sampling(::Type{SystematicSampling}, volume::AbstractVector{<:Real}, plot_area::Real, total_area::Real; kwargs...)

Performs systematic sampling for forest inventory analysis using the method of successive differences.

# Description

In systematic sampling, plots are laid out at a fixed interval along one or more lines
(transects) instead of being chosen at random, which is often far cheaper to execute in
the field. Because consecutive plots tend to be more alike than two random plots, its
precision cannot be estimated with the simple-random-sampling variance formula; the
method of successive differences below is the standard substitute (Cochran, 1977, §8.6).

# Arguments

- `volume`: vector of plot volumes, **in the order the plots were measured along the line(s)**, as a `Vol` quantity. Plain numbers are taken to be `u"m^3"`.
- `plot_area`: the area of each sample plot, as an `Area` quantity. A plain number is taken to be hectares.
- `total_area`: the total area of the forest or stand, as an `Area` quantity. A plain number is taken to be hectares.
- `line::Union{AbstractVector,Nothing}=nothing`: an optional grouping vector, one entry per plot, identifying which line/transect each plot belongs to (same order as `volume`). When given, the difference between the last plot of one line and the first plot of the next is excluded, since they are not actually adjacent on the ground. `nothing` (the default) treats every plot as a single line.
- `e::Real=10`: desired relative error margin as a percentage, only used to report whether it was met (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- `DataFrame`: a single-row report, one column per statistic (same shape as [`SimpleCasualSampling`](@ref)). Columns:
  - `vm`, `cv`, `s2m`, `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`, `pop`, `f`: as in [`SimpleCasualSampling`](@ref).
  - `ereq`: the target relative error `e`, echoed back for comparison against `relerr` (this design reports no required sample size — see Technical description below).
  - `k`: number of lines (transects).
  - `n`: number of measured plots.
  - `N`: number of possible plots in the population.
  - `area`: total area.

# Mathematical basis

For `k` lines totalling `n` plots, the variance of the mean is estimated from the sum of
squared successive differences within each line:
```math
s^2_{\\bar{x}} = \\frac{\\sum (x_i - x_{i+1})^2}{2n(n-k)}\\left(1 - \\frac{n}{N}\\right)
```
which reduces to the classic single-line formula `Σ(xᵢ-xᵢ₊₁)² / (2n(n-1))` when `k = 1`.

# Technical description

Systematic sampling has no general closed-form sample-size formula the way simple or
stratified sampling do — its precision depends on how the underlying variable happens to
be arranged along the line, not just on `n` — so, matching the source this design was
ported from, this function reports whether the achieved relative error (`relerr`) meets
the target (`ereq`) instead of a required plot count.

# Examples

```julia-repl
julia> v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1];

julia> report = sampling(SystematicSampling, v, 0.05, 10; e=10, α=0.95)
julia> report.vm
441.6909090909091 m^3
julia> report.relerr <= report.ereq
true
```
"""
struct SystematicSampling <: SamplingDesign end

"""
    MultistartSystematicSampling <: SamplingDesign

    sampling(::Type{MultistartSystematicSampling}, start::Symbol, volume::Symbol, plot_area::Area, total_area::Area,
             data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    sampling(::Type{MultistartSystematicSampling}, start::Symbol, volume::Symbol, plot_area::Real, total_area::Real,
             data::AbstractDataFrame; kwargs...)

Performs systematic sampling with multiple random starts for forest inventory analysis, generic to any number of starts.

# Description

A single systematic line (see [`SystematicSampling`](@ref)) gives no internal replication
to estimate its own precision from first principles — the successive-differences method
is an approximation. Laying out several independent systematic lines instead, each at its
own random starting point along the sampling interval, turns the design into a genuine
probability sample: every line is itself a systematic subsample of `M` plots, and the `n`
lines play the role of primary sampling units drawn from the `N` possible starting points.

# Arguments

- `start::Symbol`: name of the column identifying which random start (systematic line) each plot belongs to.
- `volume::Symbol`: name of the volume column, as a `Vol` quantity. Plain numbers are taken to be `u"m^3"`.
- `plot_area`: the area of each individual plot, as an `Area` quantity. A plain number is taken to be hectares.
- `total_area`: the total area of the forest or stand, as an `Area` quantity. A plain number is taken to be hectares.
- `data::AbstractDataFrame`: the plot-level data, one row per plot.
- `e::Real=10`: desired relative error margin as a percentage (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `startTable` (per-start descriptive statistics) and
  `resultTable` (one row, one column per statistic — same shape and column names as
  [`ClusterSampling`](@ref), since both share the same estimator).

# Mathematical basis

Statistically identical to [`ClusterSampling`](@ref)'s between/within decomposition —
each start is a primary unit of `M` plots — with the population size `N` given the
systematic-sampling interpretation of "number of possible random starts" rather than
"number of possible spatial clusters": `N = total_area / (plot_area × M)`, the number of
non-overlapping systematic subsamples that tile the population.

# Technical description

The only reason this is a separate design from [`ClusterSampling`](@ref) rather than a
mode flag is field procedure, not statistics: a "start" is a systematically-spread
subsample covering the whole area, while a "cluster" is a spatially contiguous group of
plots. Their estimators share the same variance decomposition, so both call the same
internal `_clustervariance` engine.

# Examples

```julia-repl
julia> using DataFrames

julia> data = DataFrame(
         start=repeat(1:4, inner=5),
         volume=[18.2, 19.1, 17.8, 18.9, 19.4, 22.4, 23.1, 21.9, 22.8, 23.4,
                 15.1, 16.0, 14.8, 15.6, 15.9, 20.5, 21.2, 19.8, 20.6, 21.0],
       );

julia> report = sampling(MultistartSystematicSampling, :start, :volume, 0.02, 15, data);
julia> resultTable(report).vm
19.375 m^3
```
"""
struct MultistartSystematicSampling <: SamplingDesign end

"""
    ClusterSampling <: SamplingDesign

    sampling(::Type{ClusterSampling}, cluster::Symbol, volume::Symbol, plot_area::Area, total_area::Area,
             data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    sampling(::Type{ClusterSampling}, cluster::Symbol, volume::Symbol, plot_area::Real, total_area::Real,
             data::AbstractDataFrame; kwargs...)

Performs one-stage cluster sampling for forest inventory analysis, generic to any number of clusters.

# Description

In cluster sampling the sampling unit is a cluster of `M` neighboring plots (a
"conglomerate") rather than an individual plot — cheaper to lay out in the field, at the
cost of the within-cluster homogeneity typically inflating the variance of the mean
relative to an equivalent simple random sample. Every cluster is assumed to contain the
same number of plots `M` (Cochran, 1977, §9.3, "clusters of equal size").

# Arguments

- `cluster::Symbol`: name of the column identifying which cluster each plot belongs to.
- `volume::Symbol`: name of the volume column, as a `Vol` quantity. Plain numbers are taken to be `u"m^3"`.
- `plot_area`: the area of each individual plot (secondary unit), as an `Area` quantity. A plain number is taken to be hectares.
- `total_area`: the total area of the forest or stand, as an `Area` quantity. A plain number is taken to be hectares.
- `data::AbstractDataFrame`: the plot-level data, one row per plot.
- `e::Real=10`: desired relative error margin as a percentage (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `clusterTable` (per-cluster descriptive statistics) and
  `resultTable` (one row, one column per statistic). `resultTable` columns:
  - `vm`, `cv`, `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`, `pop`, `f`: as in [`SimpleCasualSampling`](@ref), using the cluster-sampling mean/variance.
  - `s2w`, `s2b`, `s2`: within-cluster, between-cluster, and total variance (per plot).
  - `icc`: intraclass correlation coefficient `ρ`.
  - `M`: number of secondary units (plots) per cluster.
  - `n`, `nreq`, `nmiss`: measured, required, and missing number of clusters.
  - `N`: number of possible clusters in the population.
  - `area`: total area.

# Mathematical basis

With `n` sampled clusters of `M` plots each out of `N` possible clusters in the
population, the within- and between-cluster variance (per plot) are:
```math
s^2_w = \\overline{s^2_j} \\qquad
s^2_b = \\frac{M\\sum_j (\\bar{x}_j-\\bar{x})^2/(n-1) - s^2_w}{M}
```
and the variance of the overall mean, with finite-population correction `f = 1-n/N`:
```math
s^2_{\\bar{x}} = f\\frac{s^2_b}{n} + \\frac{s^2_w}{nM}
```
The intraclass correlation coefficient `ρ = s²_b/(s²_b+s²_w)` measures how similar plots
within the same cluster are: `ρ = 0` recovers the simple-random-sampling variance exactly.

# Technical description

`ClusterSampling` with `M = 1` (one plot per cluster) is numerically identical to
[`SimpleCasualSampling`](@ref) — clustering stops mattering once there is nothing left
to cluster. The population size `N` is the number of possible *clusters*, not plots:
`total_area / (plot_area × M)`.

# Examples

```julia-repl
julia> using DataFrames

julia> data = DataFrame(
         cluster=repeat(1:6, inner=4),
         volume=[18.2, 19.1, 17.8, 18.9, 22.4, 23.1, 21.9, 22.8, 15.1, 16.0, 14.8, 15.6,
                 27.3, 28.1, 26.9, 27.8, 19.8, 20.5, 19.1, 20.0, 24.5, 25.2, 23.9, 24.8],
       );

julia> report = sampling(ClusterSampling, :cluster, :volume, 0.02, 15, data);
julia> resultTable(report).vm
21.4 m^3
julia> clusterTable(report)
```
"""
struct ClusterSampling <: SamplingDesign end

"""
    HorizontalPointSampling <: SamplingDesign

    sampling(::Type{HorizontalPointSampling}, point::Symbol, diameter::Symbol, volume::Symbol, baf::Real,
             npoints::Integer, point_area::Area, total_area::Area, data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    sampling(::Type{HorizontalPointSampling}, point::Symbol, diameter::Symbol, volume::Symbol, baf::Real,
             npoints::Integer, point_area::Real, total_area::Real, data::AbstractDataFrame; kwargs...)

Performs horizontal point sampling (Bitterlich's angle-count method) for forest inventory analysis.

# Description

In horizontal point sampling, an observer stands at each sample point and uses an angle
gauge with a fixed basal area factor (`baf`) to decide, tree by tree, whether it is "in"
the count: a tree is in whenever its diameter, viewed from the point, subtends an angle
at least as wide as the critical angle the instrument defines. Every counted tree then
represents exactly `baf` (in `m²/ha`, or `ft²/ac` for imperial diameters) of basal area
per hectare, **regardless of its own size** — the defining property of the method, and
the reason it needs no fixed plot radius at all: the "plot" a tree belongs to grows with
the tree itself. This makes horizontal point sampling far faster to execute than a
fixed-area design of comparable precision, at the cost of overweighting large trees
relative to small ones — a bias only in the sense that it must be corrected for, not a
flaw, since the correction is exact.

# Arguments

- `point::Symbol`: name of the column identifying which sample point each tree was counted at.
- `diameter::Symbol`: name of the diameter column, as a `Len` quantity. Plain numbers are taken to be `u"cm"`.
- `volume::Symbol`: name of the volume column, as a `Vol` quantity. Plain numbers are taken to be `u"m^3"`.
- `baf::Real`: the instrument's basal area factor — the basal area (per hectare) that every counted tree represents, conventionally a bare number (e.g. `2.0` for `2 m²/ha` in the metric system, or the `ft²/ac` equivalent for imperial diameters).
- `npoints::Integer`: total number of sample points actually visited, **including any with zero counted trees**. Must be at least as large as the number of distinct points appearing in `data[!, point]`; the difference is treated as that many additional zero-volume observations, exactly as they must be for the sample mean to stay unbiased.
- `point_area`: the area each sample point represents — typically the grid spacing area between points on a systematic layout (not a literal fixed-radius plot, which this design has none of) — as an `Area` quantity. A plain number is taken to be hectares.
- `total_area`: the total area of the forest or stand being inventoried, as an `Area` quantity. A plain number is taken to be hectares.
- `e::Real=10`: desired relative error margin as a percentage (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `pointTable` (per-point tree count, basal area/ha, and volume/ha — only points with at least one counted tree) and `resultTable` (one row, one column per statistic):
  - `vha`: estimated volume per hectare (the primary estimator; every error/CI/sample-size statistic below is computed from it, exactly as `vm` drives them in [`SimpleCasualSampling`](@ref)).
  - `cv`, `s2m`, `se`, `abserr`, `relerr`: as in [`SimpleCasualSampling`](@ref), computed over the `npoints` per-point volume/ha values (including the zero-tree points).
  - `vtotal`, `cilower`, `ciupper`: total volume and its confidence interval.
  - `Gha`, `Gtotal`: mean basal area per hectare and its stand total — the method's most direct output, needing no individual tree volumes at all.
  - `pop`, `f`: population classification and finite-population correction factor.
  - `n`, `nreq`, `nmiss`: measured, required, and missing number of points.
  - `N`: number of possible points in the population, `total_area / point_area`.
  - `area`: total area.

# Mathematical basis

For a counted tree with basal area `g` (`ForestFoundations.basalarea(diameter)`), its individual expansion factor is
```math
EF = \\frac{BAF}{g}
```
in trees/ha — large trees have small `g` relative to `BAF` and so represent fewer
trees/ha, exactly compensating for being easier to count from farther away. Basal area
and volume per hectare at a point follow by summing each in-tree's contribution:
```math
G_{ha} = \\sum_{i} EF_i \\, g_i = n_{trees} \\cdot BAF \\qquad v_{ha} = \\sum_i EF_i \\, v_i
```
(`Gₕₐ = n_trees · BAF` holds exactly, since `EFᵢ gᵢ = BAF` for every tree by construction
— a useful sanity check, and the reason basal area alone can be estimated from tree
*counts*, without ever measuring a diameter.) The per-point `vha` values, one per sample
point (zero for points with no in-trees), are then treated as `npoints` independent
observations of a simple random sample, using the same variance/error/required-size
machinery as [`SimpleCasualSampling`](@ref) — the only difference from a fixed-area
design is that each observation is already a per-hectare rate, so the stand total is
`vha × total_area` directly rather than expanding a per-plot volume by a plot count.

# Technical description

Population size `N = total_area / point_area` and the finite-population rule follow
[`SimpleCasualSampling`](@ref) exactly, with `point_area` standing in for the area a
fixed plot would have occupied — the grid cell each point notionally controls. See
Husch, Beers, & Kershaw (2003, Ch. 15) or Avery & Burkhart (2001, Ch. 8) for the full
derivation of angle-count sampling; the method itself is due to Bitterlich (1948).

# Examples

```julia-repl
julia> using DataFrames

julia> data = DataFrame(
         point    = [1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5],
         diameter = [25.0, 30.0, 20.0, 28.0, 22.0, 35.0, 30.0, 25.0, 20.0, 22.0, 30.0, 28.0, 26.0],
         volume   = [0.35, 0.55, 0.22, 0.48, 0.28, 0.85, 0.55, 0.35, 0.22, 0.28, 0.55, 0.48, 0.40],
       );   # point 6 was visited but had no "in" trees -- not a row in `data`

julia> report = sampling(HorizontalPointSampling, :point, :diameter, :volume, 2.0, 6, 0.1, 10, data);

julia> resultTable(report).vha
1-element Vector{Quantity{Float64, 𝐋, Unitful.FreeUnits{(ha^-1, m^3), 𝐋, nothing}}}:
 32.766571189807216 m^3 ha^-1

julia> resultTable(report).Gha[1]   # from tree counts alone: 13 trees × 2 m²/ha / 6 points
4.333333333333333 m^2 ha^-1

julia> pointTable(report)
5×4 DataFrame
 Row │ point  ntrees  Gha            vha
     │ Int64  Int64   Quantity…      Quantity…
─────┼─────────────────────────────────────────────────
   1 │     1       3  6.0 m^2 ha^-1  43.8277 m^3 ha^-1
   2 │     2       2  4.0 m^2 ha^-1  30.3224 m^3 ha^-1
   3 │     3       4  8.0 m^2 ha^-1  61.4972 m^3 ha^-1
   4 │     4       1  2.0 m^2 ha^-1  14.7317 m^3 ha^-1
   5 │     5       3  6.0 m^2 ha^-1  46.2204 m^3 ha^-1
```

Point 6 was visited but had no "in" trees; it never appears in `pointTable` (only
points with at least one counted tree do), but it is still one of the `npoints=6`
observations behind `resultTable` — dropping it instead would inflate `vha`.
"""
struct HorizontalPointSampling <: SamplingDesign end

"""
    TwoStageSampling <: SamplingDesign

    sampling(::Type{TwoStageSampling}, primary::Symbol, volume::Symbol, plot_area::Area, N::Integer, M::Integer,
             data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    sampling(::Type{TwoStageSampling}, primary::Symbol, volume::Symbol, plot_area::Real, N::Integer, M::Integer,
             data::AbstractDataFrame; kwargs...)

Performs two-stage sampling for forest inventory analysis, generic to any number of primary units.

# Description

Two-stage sampling draws `n` primary units from a population of `N` (e.g. stands or
blocks), then sub-samples `m` secondary units (plots) from within each of the `M` possible
secondary units of every drawn primary — unlike [`ClusterSampling`](@ref), which
measures every secondary unit inside a drawn cluster. Sub-sampling trades some precision
for a further reduction in fieldwork, since not every drawn primary unit needs a full
census.

# Arguments

- `primary::Symbol`: name of the column identifying which primary unit each plot belongs to.
- `volume::Symbol`: name of the volume column, as a `Vol` quantity. Plain numbers are taken to be `u"m^3"`.
- `plot_area`: the area of each plot (secondary unit), as an `Area` quantity, used only to report the per-hectare volume. A plain number is taken to be hectares.
- `N::Integer`: total number of primary units in the population.
- `M::Integer`: total number of secondary units per primary unit in the population (assumed equal across primary units).
- `data::AbstractDataFrame`: the plot-level data, one row per measured plot; every sampled primary unit must contribute the same number `m` of measured plots.
- `e::Real=10`: desired relative error margin as a percentage (default 10%).
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `primaryTable` (per-primary-unit descriptive statistics)
  and `resultTable` (one row, one column per statistic). `resultTable` columns:
  - `vm`, `cv`, `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`, `pop`, `f`: as in [`SimpleCasualSampling`](@ref), using the two-stage mean/variance.
  - `s2w`, `s2b`: within-primary and between-primary variance.
  - `m`: number of secondary units measured per primary unit.
  - `M`: number of secondary units per primary unit in the population.
  - `n`, `nreq`, `nmiss`: measured, required, and missing number of primary units.
  - `N`: number of possible primary units in the population.

# Mathematical basis

With sampling fractions `f₁ = n/N` (primary stage) and `f₂ = m/M` (secondary stage):
```math
s^2_{\\bar{x}} = \\frac{1-f_1}{n}s^2_b + \\frac{f_1(1-f_2)}{nm}s^2_w
```
Setting `m = M` (every drawn primary fully censused) removes the secondary-stage
contribution entirely and this reduces to [`ClusterSampling`](@ref)'s formula.

# Technical description

The required number of primary units is solved the same way as every other design in
this module (see `_requiredsamplesize`), except the finite-population correction term
uses the *population* secondary-unit count `M` while the point estimate itself uses the
*sampled* count `m` — the only design in this module where those two differ.

# Examples

```julia-repl
julia> using DataFrames

julia> data = DataFrame(
         primary=repeat(1:5, inner=3),
         volume=[18.2, 19.1, 17.8, 22.4, 23.1, 21.9, 15.1, 16.0, 14.8, 27.3, 28.1, 26.9, 19.8, 20.5, 19.1],
       );

julia> report = sampling(TwoStageSampling, :primary, :volume, 0.02, 40, 6, data);
julia> resultTable(report).vm
20.673333333333336 m^3
```
"""
struct TwoStageSampling <: SamplingDesign end

"""
    IndependentOccasionsSampling <: SamplingDesign

    sampling(::Type{IndependentOccasionsSampling}, volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Vol}, plot_area::Area,
             N1::Integer, N2::Integer; α::Real=0.95)
    sampling(::Type{IndependentOccasionsSampling}, volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Real}, plot_area::Real,
             N1::Integer, N2::Integer; kwargs...)

Estimates growth between two inventory occasions using two independent (unmatched) samples, generic to any sample sizes.

# Description

The simplest of the successive-occasions designs: occasion 1 and occasion 2 are sampled
independently, with no attempt to remeasure the same plots. This is easy to execute but
the least efficient at detecting growth, since none of the natural plot-to-plot
correlation between the two occasions is exploited — see [`CompleteReplacementSampling`](@ref)
and [`PartialReplacementSampling`](@ref) for designs that remeasure some or all plots.

# Arguments

- `volume1`, `volume2`: independent vectors of plot volumes from occasion 1 and occasion 2 respectively, as `Vol` quantities (need not be the same length). Plain numbers are taken to be `u"m^3"`.
- `plot_area`: the area of each plot, as an `Area` quantity, used only to report the per-hectare volume. A plain number is taken to be hectares.
- `N1`, `N2`: total number of possible plots in the population at occasion 1 and occasion 2 respectively.
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `occasion1`, `occasion2` and `change`, each a single-row,
  one-column-per-statistic `DataFrame`:
  - `occasion1`/`occasion2` columns: `vm` (mean), `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`, `n`, `N`.
  - `change` columns: `gm` (mean growth), `se`, `abserr`, `relerr`, `gtotal` (total growth), `cilower`, `ciupper`.

# Mathematical basis

Since the two samples are independent, the variance of the estimated growth is simply the
sum of each occasion's variance of the mean, with no covariance term:
```math
\\bar{g} = \\bar{x}_2 - \\bar{x}_1 \\qquad s^2_{\\bar{g}} = s^2_{\\bar{x}_1} + s^2_{\\bar{x}_2}
```
using `df = (n_1-1)+(n_2-1)-1` degrees of freedom for the growth confidence interval.

# Examples

```julia-repl
julia> v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0];
julia> v2 = [24.1, 23.8, 22.9, 25.6, 24.0];

julia> report = sampling(IndependentOccasionsSampling, v1, v2, 0.05, 200, 200);
julia> change(report).gm
3.913333333333334 m^3
```
"""
struct IndependentOccasionsSampling <: SamplingDesign end

"""
    CompleteReplacementSampling <: SamplingDesign

    sampling(::Type{CompleteReplacementSampling}, volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Vol}, plot_area::Area,
             N1::Integer, N2::Integer; α::Real=0.95)
    sampling(::Type{CompleteReplacementSampling}, volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Real}, plot_area::Real,
             N1::Integer, N2::Integer; kwargs...)

Estimates growth between two inventory occasions by remeasuring the exact same plots on both occasions, generic to any sample size.

# Description

The same `n` permanent plots are measured at occasion 1 and remeasured at occasion 2
("complete/total replacement" in the sampling-on-successive-occasions terminology). Since
`volume1[i]` and `volume2[i]` refer to the same plot, the two occasions are correlated
rather than independent — usually positively, since a plot with above-average volume at
occasion 1 tends to remain above average at occasion 2 — and exploiting that correlation
tightens the growth estimate compared to [`IndependentOccasionsSampling`](@ref).

# Arguments

- `volume1`, `volume2`: matched vectors of plot volumes at occasion 1 and occasion 2 (`volume1[i]`/`volume2[i]` are the same plot), as `Vol` quantities of equal length. Plain numbers are taken to be `u"m^3"`.
- `plot_area`: the area of each plot, as an `Area` quantity, used only to report the per-hectare volume. A plain number is taken to be hectares.
- `N1`, `N2`: total number of possible plots in the population at occasion 1 and occasion 2 respectively.
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `occasion1`, `occasion2` and `change` — same column
  layout as [`IndependentOccasionsSampling`](@ref).

# Mathematical basis

The matched design subtracts twice the sample covariance between the two occasions from
the variance of the growth estimate:
```math
\\bar{g} = \\bar{x}_2 - \\bar{x}_1 \\qquad
s^2_{\\bar{g}} = s^2_{\\bar{x}_1} + s^2_{\\bar{x}_2} - \\frac{2\\,\\mathrm{Cov}(x_1,x_2)}{n}
```
using `df = 2(n-1)-1` degrees of freedom for the growth confidence interval. Positive
correlation between occasions makes `s²_ḡ` smaller than the independent-samples case.

# Examples

```julia-repl
julia> v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 24.6, 17.3];
julia> v2 = [22.1, 25.8, 23.4, 21.0, 26.6, 20.9, 26.0, 22.5];

julia> report = sampling(CompleteReplacementSampling, v1, v2, 0.05, 200, 200);
julia> change(report).gm
3.1750000000000007 m^3
```
"""
struct CompleteReplacementSampling <: SamplingDesign end

"""
    PartialReplacementSampling <: SamplingDesign

    sampling(::Type{PartialReplacementSampling}, volume1::AbstractVector{<:Union{Missing,Vol}}, volume2::AbstractVector{<:Union{Missing,Vol}},
             plot_area::Area, N::Integer; α::Real=0.95)
    sampling(::Type{PartialReplacementSampling}, volume1::AbstractVector{<:Union{Missing,Real}}, volume2::AbstractVector{<:Union{Missing,Real}},
             plot_area::Real, N::Integer; kwargs...)

Estimates growth between two inventory occasions by remeasuring only part of the plot network, generic to any sample size.

# Description

A middle ground between [`CompleteReplacementSampling`](@ref) (remeasure every plot)
and [`DoubleSampling`](@ref) (remeasure a fixed regression subset): occasion 1 has `u`
temporary and `m` permanent plots; at occasion 2 the `m` permanent plots are remeasured
and `v` brand-new temporary plots are added. Three groups of plots therefore appear across
the two vectors:

- occasion-1-only ("unmatched" — `volume1[i]` present, `volume2[i]` missing);
- measured on both occasions ("matched" — both present);
- occasion-2-only ("new temporary" — `volume1[i]` missing, `volume2[i]` present).

# Arguments

- `volume1`, `volume2`: occasion-1 and occasion-2 volumes, same length, one entry per plot, `missing` where that plot wasn't measured on that occasion (see grouping above). Plain numbers are taken to be `u"m^3"`.
- `plot_area`: the area of each plot, as an `Area` quantity, used only to report the per-hectare volume. A plain number is taken to be hectares.
- `N::Integer`: total number of possible plots in the population.
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `occasion1` (simple-random-sampling report over the `u+m`
  occasion-1 plots, same columns as [`IndependentOccasionsSampling`](@ref)'s occasion
  tables), `occasion2` (the combined estimate) and `change` (the estimated growth).
  `occasion2` columns:
  - `vm`: the combined (regression + new-temporary) occasion-2 mean estimate.
  - `vreg`: the regression-only estimate.
  - `vtemp`: the new-temporary-only mean.
  - `c`: the inverse-variance regression weight.
  - `b`: the regression slope.
  - `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`: as in [`SimpleCasualSampling`](@ref).
  - `m`, `u`, `v`: matched, unmatched, and new-temporary plot counts.
  - `N`: number of possible plots in the population.

# Mathematical basis

Two independent estimates of the occasion-2 mean are combined by inverse-variance
weighting: a regression estimate `Ŷ_reg` built the same way as in [`DoubleSampling`](@ref)
from the matched plots (using `_olsslope` for the slope), and the simple mean `Ȳᵥ` of the
new temporary plots:
```math
\\hat{Y}_{reg} = \\bar{y}_m + b(\\bar{x}_{u+m} - \\bar{x}_m) \\qquad
\\hat{Y} = c\\,\\hat{Y}_{reg} + (1-c)\\bar{Y}_v \\qquad
c^* = \\frac{s^2_{\\bar{Y}_v}}{s^2_{\\hat{Y}_{reg}} + s^2_{\\bar{Y}_v}}
```
`c*` is the variance-minimizing weight, giving `s²_Ŷ = c*·s²_Ŷreg` at the optimum.

# Technical description

Setting `v = 0` (no new temporary plots) degenerates the design toward
[`DoubleSampling`](@ref); setting `u = v = 0` (every plot matched) degenerates it toward
[`CompleteReplacementSampling`](@ref) — this design spans the space between them, and is
the most field-flexible of the four successive-occasions designs in this module.

# Examples

```julia-repl
julia> volume1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, missing, missing];
julia> volume2 = [missing, missing, 23.4, 24.0, 26.6, 22.9, 27.5, 25.2, 21.8];

julia> report = sampling(PartialReplacementSampling, volume1, volume2, 0.05, 200);
julia> occasion2(report).vm
24.34590452296151 m^3
```
"""
struct PartialReplacementSampling <: SamplingDesign end

"""
    DoubleSampling <: SamplingDesign

    sampling(::Type{DoubleSampling}, volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Union{Missing,Vol}}, plot_area::Area,
             N::Integer; α::Real=0.95)
    sampling(::Type{DoubleSampling}, volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Union{Missing,Real}}, plot_area::Real,
             N::Integer; kwargs...)

Estimates growth between two inventory occasions using a regression estimator over a large first-occasion sample and a smaller remeasured (permanent) subset, generic to any sample size.

# Description

`n` plots are measured at occasion 1; only `m ≤ n` of them ("permanent" plots) are
remeasured at occasion 2 — the rest are "temporary", occasion-1-only plots that exist
purely to pin down `x̄₁` more precisely at low cost. The occasion-2 mean is then estimated
by regressing the permanent plots' occasion-2 volume on their occasion-1 volume and
applying that regression to the full, larger occasion-1 sample (Cochran, 1977, §12.9).

# Arguments

- `volume1::AbstractVector{<:Vol}`: occasion-1 volume for every plot (`n` entries), as a `Vol` quantity. Plain numbers are taken to be `u"m^3"`.
- `volume2::AbstractVector{<:Union{Missing,Vol}}`: occasion-2 volume, one entry per plot in the same order as `volume1`; `missing` marks a temporary plot that was not remeasured.
- `plot_area`: the area of each plot, as an `Area` quantity, used only to report the per-hectare volume. A plain number is taken to be hectares.
- `N::Integer`: total number of possible plots in the population.
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `occasion1` (simple-random-sampling report over all `n`
  plots, same columns as [`IndependentOccasionsSampling`](@ref)'s occasion tables),
  `occasion2` (the regression estimate) and `change` (the estimated growth). `occasion2`
  columns:
  - `vm`: the regression-estimated occasion-2 mean.
  - `b`: the regression slope.
  - `s2yx`: residual variance of the regression.
  - `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`: as in [`SimpleCasualSampling`](@ref).
  - `m`: number of permanent (remeasured) plots.
  - `ntemp`: number of temporary plots (`n - m`).
  - `N`: number of possible plots in the population.

# Mathematical basis

The regression slope `b` is solved via the normal equations `(X'X)β = X'y` over the
permanent plots (see `_olsslope`), the same transposed-design-matrix approach used for
the ANOVA in [`StratifiedSampling`](@ref):
```math
\\hat{y}_{reg} = \\bar{y}_m + b(\\bar{x}_n - \\bar{x}_m)
```
where `x̄ₙ` is the occasion-1 mean over all `n` plots and `x̄ₘ`/`ȳₘ` are the occasion-1/2
means over only the `m` permanent plots. Its variance combines the regression's residual
variance `S²yx` (from the `m` permanent plots) with the extra precision gained from the
larger first-phase sample:
```math
s^2_{\\hat{y}_{reg}} = \\frac{S^2_{yx}}{m} + \\frac{S^2_y - S^2_{yx}}{n}
```

# Examples

```julia-repl
julia> volume1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, 17.6, 20.8, 21.9];
julia> volume2 = [22.1, 25.8, 23.4, missing, 26.6, missing, 27.5, missing, missing, 26.0];

julia> report = sampling(DoubleSampling, volume1, volume2, 0.05, 200);
julia> occasion2(report).vm
24.441857923497274 m^3
```
"""
struct DoubleSampling <: SamplingDesign end
