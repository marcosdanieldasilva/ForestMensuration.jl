"""
    taperedvolume(fit::TaperFit, dbh::Len, height::Len, hmin::Len, hmax::Len)
    taperedvolume(fit::TaperFit, dbh::Real, height::Real, hmin::Real, hmax::Real)
    taperedvolume(fit::TaperFit, dbh::Len, height::Len)
    taperedvolume(fit::TaperFit, dbh::Real, height::Real)

Volume between heights `hmin` and `hmax`, for a tree with the given
breast-height diameter `dbh` and total `height`, integrated from the fitted
stem taper curve `fit` — the curve-fitted counterpart of [`cubage`](@ref),
which integrates a *measured* section profile instead of a fitted one. The
3-argument form integrates the whole stem (`hmin=0`, `hmax=height`).

# Arguments
- `fit::TaperFit`: a stem taper model fit with `ForestModeling.fit`.
- `dbh::Len`: breast-height diameter of the tree being evaluated.
- `height::Len`: total height of the tree being evaluated.
- `hmin::Len`, `hmax::Len`: the section of the stem to integrate, `0 <= hmin < hmax <= height`.

When arguments are plain numbers, `dbh` is assumed to be in `cm` and
`height`/`hmin`/`hmax` in `m`, matching [`taperdiameter`](@ref).

# Returns
- `Vol`: the integrated volume, in `m^3`.

# Mathematical basis
```math
v = \\int_{h_{min}}^{h_{max}} \\frac{\\pi}{4} d(h)^2\\, dh
```
where `d(h)` is [`taperdiameter`](@ref)`(fit, dbh, height, h)`, integrated
numerically via `QuadGK.jl`.

# Examples
```julia-repl
julia> dbh = [20.0, 20.0, 20.0, 30.0, 30.0, 30.0]; height = [18.0, 18.0, 18.0, 22.0, 22.0, 22.0];

julia> hi = [0.3, 6.0, 14.0, 0.3, 8.0, 18.0]; di = [22.4, 15.8, 6.1, 33.6, 24.2, 8.7];

julia> ft = fit(Kozak1969(), dbh, height, hi, di);

julia> taperedvolume(ft, 25.0u"cm", 20.0u"m", 0.3u"m", 12.0u"m")
0.3170163895013209 m^3

julia> taperedvolume(ft, 25.0, 20.0)
0.4700281470314726 m^3
```
"""
function taperedvolume(fit::TaperFit, dbh::Len, height::Len, hmin::Len, hmax::Len)
  dbhval, heightval = ustrip(DUNIT, dbh), ustrip(HUNIT, height)
  hminval, hmaxval = ustrip(HUNIT, hmin), ustrip(HUNIT, hmax)
  (0 <= hminval < hmaxval <= heightval) ||
    throw(DomainError((hmin, hmax), "hmin and hmax must satisfy 0 <= hmin < hmax <= height."))
  return _taperedvolumevalue(fit, dbhval, heightval, hminval, hmaxval) * VUNIT
end

taperedvolume(fit::TaperFit, dbh::Real, height::Real, hmin::Real, hmax::Real) =
  taperedvolume(fit, dbh * DUNIT, height * HUNIT, hmin * HUNIT, hmax * HUNIT)

taperedvolume(fit::TaperFit, dbh::Len, height::Len) = taperedvolume(fit, dbh, height, zero(height), height)

taperedvolume(fit::TaperFit, dbh::Real, height::Real) =
  taperedvolume(fit, dbh * DUNIT, height * HUNIT)
