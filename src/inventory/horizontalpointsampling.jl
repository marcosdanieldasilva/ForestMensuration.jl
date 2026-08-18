# Per-point sample: sums the individual tree contributions -- basal area/ha and
# volume/ha -- for every point that had at least one "in" tree. Points with zero "in"
# trees never appear as rows in `data` and are added back as explicit zero observations
# by the caller (see `nzero` below) -- omitting them would systematically inflate the
# mean, a well-known pitfall when processing angle-count data.
function _pointtable(point::Symbol, ef::AbstractVector{<:Quantity}, g::AbstractVector{<:Area},
  vol::AbstractVector{<:Vol}, data::AbstractDataFrame)

  working = DataFrame(point=data[!, point], ef=ef, g=g, vol=vol)
  return combine(groupby(working, :point, sort=true)) do df
    (ntrees=nrow(df), Gha=sum(df.ef .* df.g), vha=sum(df.ef .* df.vol))
  end
end

"""
    horizontalpointsampling(point::Symbol, diameter::Symbol, volume::Symbol, baf::Real,
                             npoints::Integer, point_area::Area, total_area::Area,
                             data::AbstractDataFrame; e::Real=10, α::Real=0.95)
    horizontalpointsampling(point::Symbol, diameter::Symbol, volume::Symbol, baf::Real,
                             npoints::Integer, point_area::Real, total_area::Real,
                             data::AbstractDataFrame; kwargs...)

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
  - `vha`: estimated volume per hectare (the primary estimator; every error/CI/sample-size statistic below is computed from it, exactly as `vm` drives them in [`simplecasualsampling`](@ref)).
  - `cv`, `s2m`, `se`, `abserr`, `relerr`: as in [`simplecasualsampling`](@ref), computed over the `npoints` per-point volume/ha values (including the zero-tree points).
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
machinery as [`simplecasualsampling`](@ref) — the only difference from a fixed-area
design is that each observation is already a per-hectare rate, so the stand total is
`vha × total_area` directly rather than expanding a per-plot volume by a plot count.

# Technical description

Population size `N = total_area / point_area` and the finite-population rule follow
[`simplecasualsampling`](@ref) exactly, with `point_area` standing in for the area a
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

julia> report = horizontalpointsampling(:point, :diameter, :volume, 2.0, 6, 0.1, 10, data);

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
function horizontalpointsampling(point::Symbol, diameter::Symbol, volume::Symbol, baf::Real,
  npoints::Integer, point_area::Area, total_area::Area, data::AbstractDataFrame; e::Real=10, α::Real=0.95)

  baf > 0 || throw(ArgumentError("baf must be a positive value."))

  d = _asdiameter(data[!, diameter])
  vol = _asvolume(data[!, volume])
  g = basalarea.(d)
  bafunit = unit(g[1]) / referencearea(unit(d[1]))
  bafQ = baf * bafunit
  ef = bafQ ./ g

  pointtable = _pointtable(point, ef, g, vol, data)
  ncounted = nrow(pointtable)
  nzero = npoints - ncounted
  nzero >= 0 || throw(ArgumentError(
    "npoints ($npoints) is smaller than the number of distinct points with at least " *
    "one tree in data[!, point] ($ncounted)."))

  n = npoints
  n > 1 || throw(ArgumentError("At least two sampled points are required."))

  vhavalues = nzero == 0 ? pointtable.vha : vcat(pointtable.vha, fill(0.0, nzero) .* unit(pointtable.vha[1]))
  Ghavalues = nzero == 0 ? pointtable.Gha : vcat(pointtable.Gha, fill(0.0, nzero) .* unit(pointtable.Gha[1]))

  N = round(Int, ustrip(uconvert(NoUnits, total_area / point_area)))
  f = 1 - n / N
  x̅ = mean(vhavalues)
  cv = std(vhavalues) / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  if infinite
    population = "infinite"
    s²x̅ = var(vhavalues) / n
    requiredpoints = _requiredsamplesize(cv^2, n, e, α)
  else
    population = "finite"
    s²x̅ = (var(vhavalues) / n) * f
    requiredpoints = _requiredsamplesize(cv^2, n, e, α, N)
  end
  sx̅ = sqrt(s²x̅)
  missingpoints = n > requiredpoints ? 0 : requiredpoints - n

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  totalvolume = x̅ * total_area
  ciupper = totalvolume + total_area * absoluteerror
  cilower = totalvolume - total_area * absoluteerror
  Gha = mean(Ghavalues)
  Gtotal = Gha * total_area

  resulttable = DataFrame(
    vha=x̅, cv=round(cv, digits=2), s2m=s²x̅, se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vtotal=totalvolume, cilower=cilower, ciupper=ciupper,
    Gha=Gha, Gtotal=Gtotal, pop=population, f=round(f, digits=3),
    n=n, nreq=requiredpoints, nmiss=missingpoints, N=N, area=total_area,
  )

  return SamplingReport((; pointTable=pointtable, resultTable=resulttable))
end

function horizontalpointsampling(point::Symbol, diameter::Symbol, volume::Symbol, baf::Real,
  npoints::Integer, point_area::Real, total_area::Real, data::AbstractDataFrame; kwargs...)
  horizontalpointsampling(point, diameter, volume, baf, npoints, point_area * AUNIT, total_area * AUNIT, data; kwargs...)
end
