"""
    doublesampling(volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Union{Missing,Vol}}, plot_area::Area,
                   N::Integer; α::Real=0.95)
    doublesampling(volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Union{Missing,Real}}, plot_area::Real,
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
  plots, same columns as [`independentoccasionssampling`](@ref)'s occasion tables),
  `occasion2` (the regression estimate) and `change` (the estimated growth). `occasion2`
  columns:
  - `vm`: the regression-estimated occasion-2 mean.
  - `b`: the regression slope.
  - `s2yx`: residual variance of the regression.
  - `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`: as in [`simplecasualsampling`](@ref).
  - `m`: number of permanent (remeasured) plots.
  - `ntemp`: number of temporary plots (`n - m`).
  - `N`: number of possible plots in the population.

# Mathematical basis

The regression slope `b` is solved via the normal equations `(X'X)β = X'y` over the
permanent plots (see `_olsslope`), the same transposed-design-matrix approach used for
the ANOVA in [`stratifiedsampling`](@ref):
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

julia> report = doublesampling(volume1, volume2, 0.05, 200);
julia> report.occasion2.vm
24.441857923497274 m^3
```
"""
function doublesampling(volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Union{Missing,Vol}}, plot_area::Area,
  N::Integer; α::Real=0.95)

  length(volume1) == length(volume2) || throw(DimensionMismatch("volume1 and volume2 must have the same length."))
  permanentidx = findall(!ismissing, volume2)
  m = length(permanentidx)
  n = length(volume1)
  m >= 3 || throw(ArgumentError("At least three permanent (remeasured) plots are required."))

  o1 = _occasionmean(volume1, N)
  occasion1 = _occasiontable(o1, N, plot_area, α)
  vunit = unit(o1.x̅)

  x1perm = ustrip.(volume1[permanentidx])
  y2perm = ustrip.(volume2[permanentidx])
  x̄1perm, ȳ2perm = mean(x1perm), mean(y2perm)
  x̄1all = ustrip(o1.x̅)

  b = _olsslope(y2perm, x1perm)
  sstY = sum(abs2, y2perm .- ȳ2perm)
  ssxy = sum((x1perm .- x̄1perm) .* (y2perm .- ȳ2perm))
  syx² = (sstY - b * ssxy) / (m - 2)
  sy² = sstY / (m - 1)

  yreg = ȳ2perm + b * (x̄1all - x̄1perm)
  meanvar = syx² / m + (sy² - syx²) / n
  yregq, sx̅ = yreg * vunit, sqrt(meanvar) * vunit

  t = -quantile(TDist(m - 1), (1 - α) / 2)
  absoluteerror = t * sx̅
  relativeerror = absoluteerror / yregq * 100
  hectarevolume = yregq / plot_area
  totalvolume = yregq * N

  occasion2 = DataFrame(
    vm=yregq, b=round(b, digits=4), s2yx=syx² * vunit^2, se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vha=hectarevolume, vtotal=totalvolume,
    cilower=totalvolume - N * absoluteerror, ciupper=totalvolume + N * absoluteerror,
    m=m, ntemp=n - m, N=N,
  )

  growth = yregq - o1.x̅
  change = _changetable(growth, meanvar, m - 1, N, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function doublesampling(volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Union{Missing,Real}}, plot_area::Real,
  N::Integer; kwargs...)
  doublesampling(volume1 * VUNIT, _asvolume(volume2), plot_area * AUNIT, N; kwargs...)
end
