"""
    taperdiameter(fit::TaperFit, dbh::Len, height::Len, h::Len)
    taperdiameter(fit::TaperFit, dbh::Real, height::Real, h::Real)
    taperdiameter(fit::TaperFit, dbh::Len, height::Len, h::AbstractVector{<:Len})
    taperdiameter(fit::TaperFit, dbh::Real, height::Real, h::AbstractVector{<:Real})

Diameter at height `h` predicted by the fitted stem taper model `fit`, for a
tree with the given breast-height diameter `dbh` and total `height`.

# Arguments
- `fit::TaperFit`: a stem taper model fit with `ForestModeling.fit`, e.g.
  `fit(Kozak2004(), dbh, height, hi, di)` (re-exported from `ForestModeling.jl`).
- `dbh::Len`: breast-height diameter of the tree being evaluated.
- `height::Len`: total height of the tree being evaluated.
- `h`: height(s) along the stem to evaluate the diameter at.

When arguments are plain numbers, `dbh` is assumed to be in `cm` and
`height`/`h` in `m` — matching this package's default units, and the unit
`fit` is assumed to have been trained on (`ForestModeling.fit` carries no
unit convention of its own; see its docstring).

# Returns
- `Len` (or a vector of `Len`): the predicted diameter(s), in `cm`.

# Mathematical basis
```math
d(h) = dbh \\cdot d_r(h/\\text{height})
```
where `dᵣ` is `fit`'s fitted relative-diameter curve.

# Examples
```julia-repl
julia> dbh = [20.0, 20.0, 20.0, 30.0, 30.0, 30.0]; height = [18.0, 18.0, 18.0, 22.0, 22.0, 22.0];

julia> hi = [0.3, 6.0, 14.0, 0.3, 8.0, 18.0]; di = [22.4, 15.8, 6.1, 33.6, 24.2, 8.7];

julia> ft = fit(Kozak1969(), dbh, height, hi, di);

julia> taperdiameter(ft, 25.0u"cm", 20.0u"m", 7.0u"m")
19.57564935064935 cm

julia> taperdiameter(ft, 25.0, 20.0, 7.0)
19.57564935064935 cm

julia> taperdiameter(ft, 25.0, 20.0, [2.0, 7.0, 12.0])
3-element Vector{Quantity{Float64, 𝐋, Unitful.FreeUnits{(cm,), 𝐋, nothing}}}:
 22.930064935064932 cm
 19.57564935064935 cm
 15.63883116883117 cm
```
"""
function taperdiameter(fit::TaperFit, dbh::Len, height::Len, h::Len)
  d = predict(fit, ustrip(DUNIT, dbh), ustrip(HUNIT, height), ustrip(HUNIT, h))
  return d * DUNIT
end

taperdiameter(fit::TaperFit, dbh::Real, height::Real, h::Real) =
  taperdiameter(fit, dbh * DUNIT, height * HUNIT, h * HUNIT)

function taperdiameter(fit::TaperFit, dbh::Len, height::Len, h::AbstractVector{<:Len})
  d = predict(fit, ustrip(DUNIT, dbh), ustrip(HUNIT, height), ustrip.(HUNIT, h))
  return d .* DUNIT
end

taperdiameter(fit::TaperFit, dbh::Real, height::Real, h::AbstractVector{<:Real}) =
  taperdiameter(fit, dbh * DUNIT, height * HUNIT, h * HUNIT)
