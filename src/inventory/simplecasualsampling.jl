# A plot occupying a full reference area (1 ha in the metric system, 1 ac in the
# imperial one) is treated as effectively infinite regardless of the sampling fraction,
# matching the convention used by the original tool this design was ported from.
_isreferencearea(a::Area) = a == 1u"ha" || a == 1u"ac"

"""
    simplecasualsampling(volume::AbstractVector{<:Vol}, plot_area::Area, total_area::Area;
                          e::Real=10, α::Real=0.95)
    simplecasualsampling(volume::AbstractVector{<:Real}, plot_area::Real, total_area::Real; kwargs...)

Performs simple random sampling for forest inventory analysis with specified plot area and total area.

# Description

The `simplecasualsampling` function calculates various statistical parameters for a forest inventory using simple random sampling methodology. It is designed to estimate the total volume of timber within a forest area by analyzing volume measurements from sample plots.

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

julia> report = simplecasualsampling(v, 0.05, 10; e=10, α=0.95)
1×17 DataFrame
 Row │ vm          cv       s2m         se         abserr     relerr   vha             vtotal      cilower     ciupper     pop      f        n      nreq   nmiss  N      area
     │ Quantity…   Float64  Quantity…   Quantity…  Quantity…  Float64  Quantity…       Quantity…   Quantity…   Quantity…   String   Float64  Int64  Int64  Int64  Int64  Quantity…
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 441.691 m^3    10.03  168.498 m^6  12.9807 m^3  28.9227 m^3     6.55  8833.82 m^3 ha^-1  88338.2 m^3  82553.6 m^3  94122.7 m^3  finite    0.945     11      7      0    200  10 ha

julia> report.vm
441.6909090909091 m^3
```
"""
function simplecasualsampling(volume::AbstractVector{<:Vol}, plot_area::Area, total_area::Area;
  e::Real=10, α::Real=0.95)

  n = length(volume)
  n > 1 || throw(ArgumentError("At least two sampled plots are required."))

  N = round(Int, ustrip(uconvert(NoUnits, total_area / plot_area)))
  f = 1 - n / N
  x̅ = mean(volume)
  cv = std(volume) / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N) || _isreferencearea(plot_area)
  if infinite
    population = "infinite"
    s²x̅ = var(volume) / n
    requiredplots = _requiredsamplesize(cv^2, n, e, α)
  else
    population = "finite"
    s²x̅ = (var(volume) / n) * f
    requiredplots = _requiredsamplesize(cv^2, n, e, α, N)
  end
  sx̅ = sqrt(s²x̅)
  missingplots = n > requiredplots ? 0 : requiredplots - n

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevolume = x̅ / plot_area
  totalvolume = x̅ * N
  ciupper = totalvolume + N * absoluteerror
  cilower = totalvolume - N * absoluteerror

  return DataFrame(
    vm=x̅, cv=round(cv, digits=2), s2m=s²x̅, se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vha=hectarevolume, vtotal=totalvolume,
    cilower=cilower, ciupper=ciupper, pop=population, f=round(f, digits=3),
    n=n, nreq=requiredplots, nmiss=missingplots, N=N, area=total_area,
  )
end

simplecasualsampling(volume::AbstractVector{<:Real}, plot_area::Real, total_area::Real; kwargs...) =
  simplecasualsampling(volume * VUNIT, plot_area * AUNIT, total_area * AUNIT; kwargs...)
