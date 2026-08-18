# Two-stage sampling's variance decomposition, unlike one-stage cluster sampling's
# `_clustervariance`, has two independent sources of sampling error: which `n` (of `N`)
# primary units were drawn, and which `m` (of `M`) secondary units were measured within
# each of them. Cochran (1977), eq. 10.4, equal-size primaries, SRS at both stages.
function _twostagevariance(means::AbstractVector{<:Real}, vars::AbstractVector{<:Real}, m::Real, M::Real, N::Real)
  n = length(means)
  ȳ = mean(means)
  withinvar = mean(vars)
  msb = m * sum(abs2, means .- ȳ) / (n - 1)
  betweenvar = (msb - withinvar) / m
  f1, f2 = n / N, m / M
  meanvar = (1 - f1) / n * betweenvar + f1 * (1 - f2) / (n * m) * withinvar
  varterm = betweenvar + withinvar / m
  correctionvarterm = betweenvar + withinvar / M
  return (; withinvar, betweenvar, meanvar, varterm, correctionvarterm)
end

"""
    twostagesampling(primary::Symbol, volume::Symbol, plot_area::Area, N::Integer, M::Integer,
                      data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    twostagesampling(primary::Symbol, volume::Symbol, plot_area::Real, N::Integer, M::Integer,
                      data::AbstractDataFrame; kwargs...)

Performs two-stage sampling for forest inventory analysis, generic to any number of primary units.

# Description

Two-stage sampling draws `n` primary units from a population of `N` (e.g. stands or
blocks), then sub-samples `m` secondary units (plots) from within each of the `M` possible
secondary units of every drawn primary — unlike [`clustersampling`](@ref), which
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
  - `vm`, `cv`, `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`, `pop`, `f`: as in [`simplecasualsampling`](@ref), using the two-stage mean/variance.
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
contribution entirely and this reduces to [`clustersampling`](@ref)'s formula.

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

julia> report = twostagesampling(:primary, :volume, 0.02, 40, 6, data);
julia> resultTable(report).vm
20.673333333333336 m^3
```
"""
function twostagesampling(primary::Symbol, volume::Symbol, plot_area::Area, N::Integer, M::Integer,
  data::AbstractDataFrame; e::Real=10, α::Real=0.95)

  vol = _asvolume(data[!, volume])
  table = _clustertable(primary, vol, data)
  m = _equalclustersize(table.n)
  n = nrow(table)
  m <= M || throw(ArgumentError("The number of measured secondary units per primary unit (m=$m) cannot exceed the population size M=$M."))

  x̅ = mean(table.x̅)
  vunit = unit(x̅)
  v2unit = vunit^2

  variance = _twostagevariance(ustrip.(table.x̅), ustrip.(table.s²), m, M, N)
  cv = sqrt(variance.withinvar + variance.betweenvar) / ustrip(x̅) * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  E = (e / 100) * x̅
  requiredprimaries = _requiredsamplesize(variance.varterm, n, ustrip(E), α, infinite ? nothing : N;
    correctionvarterm=variance.correctionvarterm)
  missingprimaries = n > requiredprimaries ? 0 : requiredprimaries - n

  population = infinite ? "infinite" : "finite"
  sx̅ = sqrt(variance.meanvar) * vunit
  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevolume = x̅ / plot_area
  totalvolume = x̅ * N * M
  ciupper = totalvolume + N * M * absoluteerror
  cilower = totalvolume - N * M * absoluteerror

  resulttable = DataFrame(
    vm=x̅, cv=round(cv, digits=2), s2w=variance.withinvar * v2unit, s2b=variance.betweenvar * v2unit,
    se=sx̅, abserr=absoluteerror, relerr=round(relativeerror, digits=2), vha=hectarevolume,
    vtotal=totalvolume, cilower=cilower, ciupper=ciupper, pop=population, f=round(1 - n / N, digits=3),
    m=m, M=M, n=n, nreq=requiredprimaries, nmiss=missingprimaries, N=N,
  )

  return SamplingReport((; primaryTable=table, resultTable=resulttable))
end

function twostagesampling(primary::Symbol, volume::Symbol, plot_area::Real, N::Integer, M::Integer,
  data::AbstractDataFrame; kwargs...)
  twostagesampling(primary, volume, plot_area * AUNIT, N, M, data; kwargs...)
end
