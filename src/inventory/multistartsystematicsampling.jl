"""
    multistartsystematicsampling(start::Symbol, volume::Symbol, plot_area::Area, total_area::Area,
                                  data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    multistartsystematicsampling(start::Symbol, volume::Symbol, plot_area::Real, total_area::Real,
                                  data::AbstractDataFrame; kwargs...)

Performs systematic sampling with multiple random starts for forest inventory analysis, generic to any number of starts.

# Description

A single systematic line (see [`systematicsampling`](@ref)) gives no internal replication
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
  [`clustersampling`](@ref), since both share the same estimator).

# Mathematical basis

Statistically identical to [`clustersampling`](@ref)'s between/within decomposition —
each start is a primary unit of `M` plots — with the population size `N` given the
systematic-sampling interpretation of "number of possible random starts" rather than
"number of possible spatial clusters": `N = total_area / (plot_area × M)`, the number of
non-overlapping systematic subsamples that tile the population.

# Technical description

The only reason this is a separate function from `clustersampling` rather than a mode
flag is field procedure, not statistics: a "start" is a systematically-spread subsample
covering the whole area, while a "cluster" is a spatially contiguous group of plots. Their
estimators share the same variance decomposition, so both call the same internal
`_clustervariance` engine.

# Examples

```julia-repl
julia> using DataFrames

julia> data = DataFrame(
         start=repeat(1:4, inner=5),
         volume=[18.2, 19.1, 17.8, 18.9, 19.4, 22.4, 23.1, 21.9, 22.8, 23.4,
                 15.1, 16.0, 14.8, 15.6, 15.9, 20.5, 21.2, 19.8, 20.6, 21.0],
       );

julia> report = multistartsystematicsampling(:start, :volume, 0.02, 15, data);
julia> resultTable(report).vm
19.375 m^3
```
"""
function multistartsystematicsampling(start::Symbol, volume::Symbol, plot_area::Area, total_area::Area,
  data::AbstractDataFrame; e::Real=10, α::Real=0.95)

  vol = _asvolume(data[!, volume])
  table = _clustertable(start, vol, data)
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

  return SamplingReport((; startTable=table, resultTable=resulttable))
end

function multistartsystematicsampling(start::Symbol, volume::Symbol, plot_area::Real, total_area::Real,
  data::AbstractDataFrame; kwargs...)
  multistartsystematicsampling(start, volume, plot_area * AUNIT, total_area * AUNIT, data; kwargs...)
end
