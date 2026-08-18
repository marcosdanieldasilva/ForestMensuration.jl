# ==============================================================================
# TAPER APPLICATION — SHARED INTERNALS
# ==============================================================================
# Every function here works in plain Float64, in the unit `fit` was trained on
# (cm for diameters, m for heights — the package-wide default convention, see
# `units.jl`), so the public, `Unitful`-aware functions in `taperdiameter.jl`/
# `taperheight.jl`/`taperedvolume.jl`/`logassortment.jl` only need to strip
# units once at the boundary and tag the result once on the way out, instead
# of repeating that conversion inside a numerical inner loop (root-finding,
# quadrature, or the greedy bucking scan in `logassortment`).

"""
    _taperheightvalue(fit::TaperFit, dbh::Real, ht::Real, d::Real) -> Float64

Height at which the fitted taper curve narrows to diameter `d`, found by
bracketed root-finding (`Roots.jl`) over `(0, ht)`. Throws `DomainError` if
`d` is not attainable anywhere on the stem (larger than the base diameter, or
smaller than the tip diameter) — the strict validation appropriate for a
direct user query; see `_hlimitfor` for the non-throwing variant `logassortment`
uses internally.
"""
function _taperheightvalue(fit::TaperFit, dbh::Real, ht::Real, d::Real)
  lo, hi = 1e-6, ht - 1e-6
  dlo, dhi = predict(fit, dbh, ht, lo), predict(fit, dbh, ht, hi)
  f(h) = predict(fit, dbh, ht, h) - d
  sign(dlo - d) == sign(dhi - d) && throw(DomainError(d,
    "diameter $d is outside the range this taper model produces for a tree of dbh=$dbh, " *
    "height=$ht (base: $dlo, tip: $dhi)."))
  return Roots.find_zero(f, (lo, hi), Roots.Bisection())
end

"""
    _taperedvolumevalue(fit::TaperFit, dbh::Real, ht::Real, hmin::Real, hmax::Real) -> Float64

Volume (in `m³`, since `dbh`/heights are cm/m) between `hmin` and `hmax`,
integrating the taper curve's cross-sectional area (`π·d(h)²/40000`, the
standard cm-diameter-to-m²-area disk-method constant) via `QuadGK.jl`.
"""
function _taperedvolumevalue(fit::TaperFit, dbh::Real, ht::Real, hmin::Real, hmax::Real)
  integrand(h) = π * predict(fit, dbh, ht, h)^2 / 40000
  v, _ = QuadGK.quadgk(integrand, hmin, hmax)
  return v
end

"""
    _hlimitfor(fit::TaperFit, dbh::Real, ht::Real, targetd::Real) -> Float64

Non-throwing counterpart to `_taperheightvalue`, used by `logassortment`'s
greedy bucking loop: the height at which the stem narrows to `targetd`, or
`ht` if the stem never narrows that far (the whole remaining length still
qualifies), or (numerically) the stem base if `targetd` is unattainable
anywhere (so the calling loop's `while` condition simply never fires — no
volume/logs assigned to that product, with no exception needed to signal it).
"""
function _hlimitfor(fit::TaperFit, dbh::Real, ht::Real, targetd::Real)
  lo, hi = 1e-6, ht - 1e-6
  predict(fit, dbh, ht, lo) < targetd && return lo
  predict(fit, dbh, ht, hi) >= targetd && return ht
  f(h) = predict(fit, dbh, ht, h) - targetd
  return Roots.find_zero(f, (lo, hi), Roots.Bisection())
end
