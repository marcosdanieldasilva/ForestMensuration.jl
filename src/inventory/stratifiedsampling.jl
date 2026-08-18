# Per-stratum descriptive statistics: sample size, mean, variance, standard deviation,
# and the stratum's relative weight `p` (its share of the total population area) together
# with `p·s` and `p·s²`, the building blocks of every stratified-sampling formula below.
# Generic to any number of strata via `groupby`.
function _auxiliarytable(stratum::Symbol, vol::AbstractVector{<:Vol}, strata_area::AbstractVector{<:Area}, data::AbstractDataFrame)
  data = copy(data)
  data[!, :__volume__] = vol
  table = combine(groupby(data, stratum, sort=true)) do df
    (n=nrow(df), x̅=mean(df.__volume__), s²=var(df.__volume__), s=std(df.__volume__))
  end
  p = strata_area ./ sum(strata_area)
  table.p = p
  table.ps = p .* table.s
  table.ps² = p .* table.s²
  return table
end

# One-way ANOVA testing whether the strata means differ, via the standard sum-of-squares
# decomposition -- mathematically identical to fitting `volume ~ 1 + stratum` and reading
# off its F-test, without needing a modeling/design-matrix dependency for what is, here,
# a single categorical predictor.
function _anova(vol::AbstractVector{<:Vol}, table::AbstractDataFrame)
  n, k = length(vol), nrow(table)
  ȳ = mean(vol)
  sst = sum(abs2, vol .- ȳ)
  ssb = sum(table.n .* abs2.(table.x̅ .- ȳ))
  ssw = sst - ssb
  dfb, dfw = k - 1, n - k
  msb, msw = ssb / dfb, ssw / dfw
  # a single stratum has no between-strata comparison to make: the F-test is undefined
  # (0 numerator degrees of freedom), so the between-strata row is left blank instead.
  fstat = dfb > 0 ? msb / msw : missing
  pval = dfb > 0 ? ccdf(FDist(dfb, dfw), fstat) : missing
  DataFrame(
    "Source" => ["Between strata", "Within strata", "Total"],
    "DOF" => [dfb, dfw, n - 1],
    "SS" => [ssb, ssw, sst],
    "MS" => [dfb > 0 ? msb : missing, msw, missing],
    "F" => [fstat, missing, missing],
    "p(>F)" => [pval, missing, missing],
  )
end

# Elementwise max(required - measured, 0): how many additional plots each stratum still
# needs, never negative when a stratum already exceeds its required allocation.
_missingplots(measured::AbstractVector{<:Integer}, required::AbstractVector{<:Integer}) = max.(required .- measured, 0)

"""
    stratifiedsampling(stratum::Symbol, volume::Symbol, plot_area::Area, strata_area::AbstractVector{<:Area},
                        data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    stratifiedsampling(stratum::Symbol, volume::Symbol, plot_area::Real, strata_area::AbstractVector{<:Real},
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
  statistic — same shape as [`simplecasualsampling`](@ref)'s return value). `resultTable`
  columns:
  - `vm`, `cv`, `s2m`, `se`, `abserr`, `relerr`, `vtotal`, `cilower`, `ciupper`, `pop`, `f`: as in [`simplecasualsampling`](@ref), using the stratified mean/variance.
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
proportion to `pₕsₕ`, solved iteratively the same way as [`simplecasualsampling`](@ref)
(see `_requiredsamplesize`). The confidence interval uses a Satterthwaite-approximated
effective degrees of freedom that accounts for the strata's differing sizes and variances
(Cochran, 1977, eq. 5.35):
```math
df = \\frac{\\left(\\sum_h g_h s_h^2\\right)^2}{\\sum_h \\dfrac{g_h^2 s_h^4}{n_h - 1}}
\\qquad \\text{where} \\qquad g_h = \\frac{N_h(N_h - n_h)}{n_h}
```

# Technical description

With a single stratum this reduces exactly to [`simplecasualsampling`](@ref) — the
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

julia> report = stratifiedsampling(:stratum, :volume, 0.1, [12.0, 8.0, 20.0], data);

julia> resultTable(report)
julia> auxiliaryTable(report)
julia> anova(report)
```
"""
function stratifiedsampling(stratum::Symbol, volume::Symbol, plot_area::Area, strata_area::AbstractVector{<:Area},
  data::AbstractDataFrame; e::Real=10, α::Real=0.95)

  length(strata_area) == length(unique(data[!, stratum])) ||
    throw(ArgumentError("strata_area must have one entry per unique stratum, in sorted order."))

  vol = _asvolume(data[!, volume])
  table = _auxiliarytable(stratum, vol, strata_area, data)
  anova = _anova(vol, table)

  N = round(Int, ustrip(uconvert(NoUnits, sum(strata_area) / plot_area)))
  x̅st = sum(table.p .* table.x̅)
  s²x̅st = sum(table.ps .^ 2 ./ table.n) - sum(table.ps²) / N
  sx̅st = sqrt(s²x̅st)
  cv = sx̅st / x̅st * 100
  E = (e / 100) * x̅st
  measured = table.n
  n = sum(measured)

  infinite = _isinfinitepopulation(n, N)
  requiredtotal = _requiredsamplesize(ustrip(sum(table.ps²)), n, ustrip(E), α, infinite ? nothing : N)
  required = round.(Int, requiredtotal .* table.ps ./ sum(table.ps))
  missingplots = _missingplots(measured, required)

  Nh = ustrip.(uconvert.(NoUnits, strata_area ./ plot_area))
  gh = Nh .* (Nh .- measured) ./ measured
  dfeff = sum(gh .* table.s²)^2 / sum((gh .^ 2 .* table.s² .^ 2) ./ (measured .- 1))
  Ttab = -quantile(TDist(dfeff), (1 - α) / 2)

  population = infinite ? "infinite" : "finite"
  f = 1 - n / N
  absoluteerror = Ttab * sx̅st
  relativeerror = absoluteerror / x̅st * 100
  totalvolume = x̅st * N
  ciupper = totalvolume + N * absoluteerror
  cilower = totalvolume - N * absoluteerror

  resulttable = DataFrame(
    vm=x̅st, cv=round(cv, digits=2), s2m=s²x̅st, se=sx̅st, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vtotal=totalvolume, cilower=cilower, ciupper=ciupper,
    pop=population, f=round(f, digits=3), nh=Tuple(measured), nreqh=Tuple(required),
    nmissh=Tuple(missingplots), n=n, nreq=requiredtotal, N=N, areah=Tuple(ustrip.(strata_area)),
  )

  return SamplingReport((; anova, auxiliaryTable=table, resultTable=resulttable))
end

function stratifiedsampling(stratum::Symbol, volume::Symbol, plot_area::Real, strata_area::AbstractVector{<:Real},
  data::AbstractDataFrame; kwargs...)
  stratifiedsampling(stratum, volume, plot_area * AUNIT, strata_area * AUNIT, data; kwargs...)
end
