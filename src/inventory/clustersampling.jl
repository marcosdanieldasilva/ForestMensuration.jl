# Per-cluster sample size, mean and internal variance, one row per cluster, generic to
# any number of clusters via `groupby`.
function _clustertable(cluster::Symbol, vol::AbstractVector{<:Vol}, data::AbstractDataFrame)
  data = copy(data)
  data[!, :__volume__] = vol
  return combine(groupby(data, cluster, sort=true)) do df
    (n=nrow(df), x̅=mean(df.__volume__), s²=var(df.__volume__))
  end
end

# Validates that every cluster has the same number of secondary units and returns it;
# `_clustervariance`'s decomposition assumes equal-size primary units.
function _equalclustersize(n::AbstractVector{<:Integer})
  sizes = unique(n)
  length(sizes) == 1 ||
    throw(ArgumentError("Every cluster must contain the same number of secondary units; found sizes $(sort(sizes))."))
  return only(sizes)
end

"""
    clustersampling(cluster::Symbol, volume::Symbol, plot_area::Area, total_area::Area,
                     data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    clustersampling(cluster::Symbol, volume::Symbol, plot_area::Real, total_area::Real,
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

- [`SamplingReport`](@ref) with `cluster_table` (per-cluster descriptive statistics) and
  `result_table` (one row, one column per statistic). `result_table` columns:
  - `vm`, `cv`, `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`, `pop`, `f`: as in [`simplecasualsampling`](@ref), using the cluster-sampling mean/variance.
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

`clustersampling` with `M = 1` (one plot per cluster) is numerically identical to
[`simplecasualsampling`](@ref) — clustering stops mattering once there is nothing left
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

julia> report = clustersampling(:cluster, :volume, 0.02, 15, data);
julia> report.result_table.vm
21.4 m^3
julia> report.cluster_table
```
"""
function clustersampling(cluster::Symbol, volume::Symbol, plot_area::Area, total_area::Area,
  data::AbstractDataFrame; e::Real=10, α::Real=0.95)

  vol = _asvolume(data[!, volume])
  table = _clustertable(cluster, vol, data)
  M = _equalclustersize(table.n)
  n = nrow(table)

  N = round(Int, ustrip(uconvert(NoUnits, total_area / (plot_area * M))))
  x̅ = mean(table.x̅)
  vunit = unit(x̅)
  v2unit = vunit^2

  variance = _clustervariance(ustrip.(table.x̅), ustrip.(table.s²), M, N)
  cv = sqrt(variance.totalvar) / ustrip(x̅) * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  E = (e / 100) * x̅
  requiredplots = _requiredsamplesize(variance.varterm, n, ustrip(E), α, infinite ? nothing : N)
  missingplots = n > requiredplots ? 0 : requiredplots - n

  population = infinite ? "infinite" : "finite"
  sx̅ = sqrt(variance.meanvar) * vunit
  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevolume = x̅ / plot_area
  totalvolume = x̅ * M * N
  ciupper = totalvolume + M * N * absoluteerror
  cilower = totalvolume - M * N * absoluteerror

  resulttable = DataFrame(
    vm=x̅, cv=round(cv, digits=2), s2w=variance.withinvar * v2unit, s2b=variance.betweenvar * v2unit,
    s2=variance.totalvar * v2unit, icc=round(variance.icc, digits=4), se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vha=hectarevolume, vtotal=totalvolume, cilower=cilower,
    ciupper=ciupper, pop=population, f=round(1 - n / N, digits=3), M=M, n=n, nreq=requiredplots,
    nmiss=missingplots, N=N, area=total_area,
  )

  return SamplingReport((; cluster_table=table, result_table=resulttable))
end

function clustersampling(cluster::Symbol, volume::Symbol, plot_area::Real, total_area::Real,
  data::AbstractDataFrame; kwargs...)
  clustersampling(cluster, volume, plot_area * AUNIT, total_area * AUNIT, data; kwargs...)
end
