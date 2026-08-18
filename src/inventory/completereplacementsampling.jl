"""
    completereplacementsampling(volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Vol}, plot_area::Area,
                                 N1::Integer, N2::Integer; α::Real=0.95)
    completereplacementsampling(volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Real}, plot_area::Real,
                                 N1::Integer, N2::Integer; kwargs...)

Estimates growth between two inventory occasions by remeasuring the exact same plots on both occasions, generic to any sample size.

# Description

The same `n` permanent plots are measured at occasion 1 and remeasured at occasion 2
("complete/total replacement" in the sampling-on-successive-occasions terminology). Since
`volume1[i]` and `volume2[i]` refer to the same plot, the two occasions are correlated
rather than independent — usually positively, since a plot with above-average volume at
occasion 1 tends to remain above average at occasion 2 — and exploiting that correlation
tightens the growth estimate compared to [`independentoccasionssampling`](@ref).

# Arguments

- `volume1`, `volume2`: matched vectors of plot volumes at occasion 1 and occasion 2 (`volume1[i]`/`volume2[i]` are the same plot), as `Vol` quantities of equal length. Plain numbers are taken to be `u"m^3"`.
- `plot_area`: the area of each plot, as an `Area` quantity, used only to report the per-hectare volume. A plain number is taken to be hectares.
- `N1`, `N2`: total number of possible plots in the population at occasion 1 and occasion 2 respectively.
- `α::Real=0.95`: confidence level (default 95%).

# Returns

- [`SamplingReport`](@ref) with `occasion1`, `occasion2` and `change` — same column
  layout as [`independentoccasionssampling`](@ref).

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

julia> report = completereplacementsampling(v1, v2, 0.05, 200, 200);
julia> change(report).gm
3.1750000000000007 m^3
```
"""
function completereplacementsampling(volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Vol}, plot_area::Area,
  N1::Integer, N2::Integer; α::Real=0.95)

  length(volume1) == length(volume2) ||
    throw(DimensionMismatch("volume1 and volume2 must have the same length: the same plots are remeasured on both occasions."))

  o1 = _occasionmean(volume1, N1)
  o2 = _occasionmean(volume2, N2)
  occasion1 = _occasiontable(o1, N1, plot_area, α)
  occasion2 = _occasiontable(o2, N2, plot_area, α)

  n = o1.n
  growth = o2.x̅ - o1.x̅
  growthcovariance = ustrip(cov(volume1, volume2))
  # with a small sample and strong occasion-to-occasion correlation, the covariance
  # estimate can occasionally overshoot enough to make this formula dip below zero;
  # clamped at zero, since a variance estimate cannot legitimately be negative.
  growthvar = max(ustrip(o1.meanvar) + ustrip(o2.meanvar) - 2 * growthcovariance / n, 0.0)
  df = 2 * (n - 1) - 1
  change = _changetable(growth, growthvar, df, N2, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function completereplacementsampling(volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Real}, plot_area::Real,
  N1::Integer, N2::Integer; kwargs...)
  completereplacementsampling(volume1 * VUNIT, volume2 * VUNIT, plot_area * AUNIT, N1, N2; kwargs...)
end
