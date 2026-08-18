"""
    taperheight(fit::TaperFit, dbh::Len, height::Len, d::Len)
    taperheight(fit::TaperFit, dbh::Real, height::Real, d::Real)

Height at which the fitted stem taper model `fit` narrows to diameter `d`,
for a tree with the given breast-height diameter `dbh` and total `height` —
the inverse of [`taperdiameter`](@ref).

# Arguments
- `fit::TaperFit`: a stem taper model fit with `ForestModeling.fit`.
- `dbh::Len`: breast-height diameter of the tree being evaluated.
- `height::Len`: total height of the tree being evaluated.
- `d::Len`: target diameter (e.g. a minimum merchantable top diameter).

When arguments are plain numbers, `dbh`/`d` are assumed to be in `cm` and
`height` in `m`, matching [`taperdiameter`](@ref).

# Returns
- `Len`: the height, in `m`, at which the taper curve equals `d`.

Throws `DomainError` if `d` is not attainable anywhere on the stem (larger
than the predicted base diameter, or smaller than the predicted tip diameter).

# Technical description
Solved by bracketed root-finding (`Roots.jl`, bisection) over `(0, height)` —
the taper curve is assumed monotonically decreasing in `h`, as every model in
this package's catalog is by construction.

# Examples
```julia-repl
julia> dbh = [20.0, 20.0, 20.0, 30.0, 30.0, 30.0]; height = [18.0, 18.0, 18.0, 22.0, 22.0, 22.0];

julia> hi = [0.3, 6.0, 14.0, 0.3, 8.0, 18.0]; di = [22.4, 15.8, 6.1, 33.6, 24.2, 8.7];

julia> ft = fit(Kozak1969(), dbh, height, hi, di);

julia> taperheight(ft, 25.0u"cm", 20.0u"m", 15.0u"cm")
12.756673511293635 m

julia> taperheight(ft, 25.0, 20.0, 15.0)
12.756673511293635 m
```
"""
function taperheight(fit::TaperFit, dbh::Len, height::Len, d::Len)
  h = _taperheightvalue(fit, ustrip(DUNIT, dbh), ustrip(HUNIT, height), ustrip(DUNIT, d))
  return h * HUNIT
end

taperheight(fit::TaperFit, dbh::Real, height::Real, d::Real) =
  taperheight(fit, dbh * DUNIT, height * HUNIT, d * DUNIT)
