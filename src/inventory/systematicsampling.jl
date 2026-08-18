# Successive differences between consecutive plots along a systematic line, skipping the
# pair that straddles two different lines when a `line` grouping is given -- generic to
# any number of lines (transects), reducing to a single running difference when there is
# only one.
function _successivedifferences(vol::AbstractVector{<:Vol}, line::Union{AbstractVector,Nothing})
  n = length(vol)
  keep = isnothing(line) ? trues(n - 1) : [line[i] == line[i+1] for i in 1:n-1]
  return [vol[i] - vol[i+1] for i in 1:n-1 if keep[i]]
end

"""
    systematicsampling(volume::AbstractVector{<:Vol}, plot_area::Area, total_area::Area;
                        line::Union{AbstractVector,Nothing}=nothing, e::Real=10, α::Real=0.95)
    systematicsampling(volume::AbstractVector{<:Real}, plot_area::Real, total_area::Real; kwargs...)

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

- `DataFrame`: a single-row report, one column per statistic (same shape as [`simplecasualsampling`](@ref)). Columns:
  - `vm`, `cv`, `s2m`, `se`, `abserr`, `relerr`, `vha`, `vtotal`, `cilower`, `ciupper`, `pop`, `f`: as in [`simplecasualsampling`](@ref).
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

julia> report = systematicsampling(v, 0.05, 10; e=10, α=0.95)
julia> report.vm
441.6909090909091 m^3
julia> report.relerr <= report.ereq
true
```
"""
function systematicsampling(volume::AbstractVector{<:Vol}, plot_area::Area, total_area::Area;
  line::Union{AbstractVector,Nothing}=nothing, e::Real=10, α::Real=0.95)

  n = length(volume)
  n > 1 || throw(ArgumentError("At least two sampled plots are required."))
  isnothing(line) || length(line) == n || throw(DimensionMismatch("line must have the same length as volume."))

  N = round(Int, ustrip(uconvert(NoUnits, total_area / plot_area)))
  k = isnothing(line) ? 1 : length(unique(line))
  x̅ = mean(volume)
  diffs = _successivedifferences(volume, line)

  infinite = _isinfinitepopulation(n, N)
  f = infinite ? 1.0 : 1 - n / N
  population = infinite ? "infinite" : "finite"
  s²x̅ = sum(abs2, diffs) / (2n * (n - k)) * f
  sx̅ = sqrt(s²x̅)
  cv = sx̅ / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevolume = x̅ / plot_area
  totalvolume = x̅ * N
  ciupper = totalvolume + N * absoluteerror
  cilower = totalvolume - N * absoluteerror

  return DataFrame(
    vm=x̅, cv=round(cv, digits=2), s2m=s²x̅, se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), ereq=e, vha=hectarevolume, vtotal=totalvolume,
    cilower=cilower, ciupper=ciupper, pop=population, f=round(f, digits=3),
    k=k, n=n, N=N, area=total_area,
  )
end

function systematicsampling(volume::AbstractVector{<:Real}, plot_area::Real, total_area::Real; kwargs...)
  systematicsampling(volume * VUNIT, plot_area * AUNIT, total_area * AUNIT; kwargs...)
end
