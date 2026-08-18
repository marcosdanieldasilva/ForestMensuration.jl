"""
    partialreplacementsampling(volume1::AbstractVector{<:Union{Missing,Vol}}, volume2::AbstractVector{<:Union{Missing,Vol}},
                                plot_area::Area, N::Integer; α::Real=0.95)
    partialreplacementsampling(volume1::AbstractVector{<:Union{Missing,Real}}, volume2::AbstractVector{<:Union{Missing,Real}},
                                plot_area::Real, N::Integer; kwargs...)

Estimates growth between two inventory occasions by remeasuring only part of the plot network, generic to any sample size.

# Description

A middle ground between [`completereplacementsampling`](@ref) (remeasure every plot)
and [`doublesampling`](@ref) (remeasure a fixed regression subset): occasion 1 has `u`
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
  occasion-1 plots, same columns as [`independentoccasionssampling`](@ref)'s occasion
  tables), `occasion2` (the combined estimate) and `change` (the estimated growth).
  `occasion2` columns:
  - `vm`: the combined (regression + new-temporary) occasion-2 mean estimate.
  - `vreg`: the regression-only estimate.
  - `vtemp`: the new-temporary-only mean.
  - `c`: the inverse-variance regression weight.
  - `b`: the regression slope.
  - `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`: as in [`simplecasualsampling`](@ref).
  - `m`, `u`, `v`: matched, unmatched, and new-temporary plot counts.
  - `N`: number of possible plots in the population.

# Mathematical basis

Two independent estimates of the occasion-2 mean are combined by inverse-variance
weighting: a regression estimate `Ŷ_reg` built the same way as in [`doublesampling`](@ref)
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
[`doublesampling`](@ref); setting `u = v = 0` (every plot matched) degenerates it toward
[`completereplacementsampling`](@ref) — this design spans the space between them, and is
the most field-flexible of the four successive-occasions designs in this module.

# Examples

```julia-repl
julia> volume1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, missing, missing];
julia> volume2 = [missing, missing, 23.4, 24.0, 26.6, 22.9, 27.5, 25.2, 21.8];

julia> report = partialreplacementsampling(volume1, volume2, 0.05, 200);
julia> occasion2(report).vm
24.34590452296151 m^3
```
"""
function partialreplacementsampling(volume1::AbstractVector{<:Union{Missing,Vol}}, volume2::AbstractVector{<:Union{Missing,Vol}},
  plot_area::Area, N::Integer; α::Real=0.95)

  length(volume1) == length(volume2) || throw(DimensionMismatch("volume1 and volume2 must have the same length."))

  idx = eachindex(volume1)
  matched = findall(i -> !ismissing(volume1[i]) && !ismissing(volume2[i]), idx)
  unmatched = findall(i -> !ismissing(volume1[i]) && ismissing(volume2[i]), idx)
  newtemp = findall(i -> ismissing(volume1[i]) && !ismissing(volume2[i]), idx)
  m, u, v = length(matched), length(unmatched), length(newtemp)
  m >= 3 || throw(ArgumentError("At least three matched (remeasured on both occasions) plots are required."))
  v >= 2 || throw(ArgumentError("At least two new temporary (occasion-2-only) plots are required."))

  o1 = _occasionmean(collect(skipmissing(volume1)), N)
  occasion1 = _occasiontable(o1, N, plot_area, α)
  vunit = unit(o1.x̅)

  x1m, y2m = ustrip.(volume1[matched]), ustrip.(volume2[matched])
  x1u = ustrip.(volume1[unmatched])
  y2v = ustrip.(volume2[newtemp])
  n1 = u + m
  x̅u = isempty(x1u) ? zero(x1m[1]) : mean(x1u)
  x̅m, ȳm = mean(x1m), mean(y2m)
  x̅combined = (u * x̅u + m * x̅m) / n1
  Sym² = var(y2m)

  b = _olsslope(y2m, x1m)
  sstY = sum(abs2, y2m .- ȳm)
  ssxy = sum((x1m .- x̅m) .* (y2m .- ȳm))
  syx² = (sstY - b * ssxy) / (m - 2)
  yreg = ȳm + b * (x̅combined - x̅m)
  varyreg = syx² / m + (Sym² - syx²) / n1

  ȳv = mean(y2v)
  varyv = var(y2v) / v * (1 - v / N)

  c = varyv / (varyreg + varyv)
  yspr = c * yreg + (1 - c) * ȳv
  varyspr = c^2 * varyreg + (1 - c)^2 * varyv

  ysprq, sx̅ = yspr * vunit, sqrt(varyspr) * vunit
  t = -quantile(TDist(m + v - 2), (1 - α) / 2)
  absoluteerror = t * sx̅
  relativeerror = absoluteerror / ysprq * 100
  hectarevolume = ysprq / plot_area
  totalvolume = ysprq * N

  occasion2 = DataFrame(
    vm=ysprq, vreg=yreg * vunit, vtemp=ȳv * vunit, c=round(c, digits=4), b=round(b, digits=4),
    se=sx̅, abserr=absoluteerror, relerr=round(relativeerror, digits=2), vha=hectarevolume,
    vtotal=totalvolume, cilower=totalvolume - N * absoluteerror, ciupper=totalvolume + N * absoluteerror,
    m=m, u=u, v=v, N=N,
  )

  growth = ysprq - x̅combined * vunit
  change = _changetable(growth, varyspr, m + v - 2, N, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function partialreplacementsampling(volume1::AbstractVector{<:Union{Missing,Real}}, volume2::AbstractVector{<:Union{Missing,Real}},
  plot_area::Real, N::Integer; kwargs...)
  partialreplacementsampling(_asvolume(volume1), _asvolume(volume2), plot_area * AUNIT, N; kwargs...)
end
