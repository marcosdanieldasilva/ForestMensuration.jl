# Greedy priority-order log bucking, generalized to any number of products.
# Walks `sed`/`minlen`/`maxlen`/`kerf` in row order (row 1 = most valuable
# product), cutting the longest permitted log for the current product until
# the stem narrows below its minimum small-end diameter (`_hlimitfor`), then
# advancing to the next product for the remainder of the stem — `h0` carries
# forward across products (never resets), so the stem is consumed strictly
# bottom-up, exactly once. This single, model-agnostic loop (parametrized
# only by `predict`/`_taperedvolumevalue`/`_hlimitfor`) is what structurally
# rules out `timbeR`'s two known bugs: there is only one `kerf` read-path (no
# risk of reading the wrong assortment-table column), and diameter/height/
# volume all route through the same `TaperFit`, so no formula can silently
# disagree with another about e.g. a Kozak2004 exponent.
function _cutlogs(fit::TaperFit, dbh::Real, ht::Real, stumpheight::Real,
  sed::AbstractVector{<:Real}, minlen::AbstractVector{<:Real}, maxlen::AbstractVector{<:Real}, kerf::AbstractVector{<:Real})

  n = length(sed)
  volumes = zeros(n)
  logcounts = zeros(Int, n)
  h0 = stumpheight
  for i in 1:n
    hlimit = _hlimitfor(fit, dbh, ht, sed[i])
    while h0 + minlen[i] <= hlimit
      seglen = min(maxlen[i], hlimit - h0)
      volumes[i] += _taperedvolumevalue(fit, dbh, ht, h0, h0 + seglen)
      logcounts[i] += 1
      h0 += seglen + kerf[i]
    end
  end
  return volumes, logcounts
end

"""
    logassortment(fit::TaperFit, dbh::Len, height::Len, products::AbstractDataFrame; stumpheight::Len=0.0u"m")
    logassortment(fit::TaperFit, dbh::Real, height::Real, products::AbstractDataFrame; stumpheight::Real=0.0)

Simulates cutting a tree's stem into logs ("sortimentos") across `N` product
classes, using the fitted stem taper model `fit` — the generalized,
bug-fixed equivalent of `timbeR`'s per-model `*_logs` functions.

# Arguments
- `fit::TaperFit`: a stem taper model fit with `ForestModeling.fit`.
- `dbh::Len`: breast-height diameter of the tree being bucked.
- `height::Len`: total height of the tree being bucked.
- `products::AbstractDataFrame`: one row per product class, in **cutting
  priority order** (row 1 = most valuable, cut first), with columns:
  - `name`: product label (any `AbstractString`-eltype column).
  - `sed`: minimum small-end (top) diameter accepted for this product.
  - `minlength`: shortest usable log length for this product.
  - `maxlength`: longest log cut per pass for this product.
  - `kerf`: stem length lost to each cut (saw kerf) after a log of this product.
  `sed` may be `Len` or a plain number (`cm`); `minlength`/`maxlength`/`kerf`
  may be `Len` or a plain number (`m`).
- `stumpheight::Len=0`: height of the stump left uncut at the base.

When `dbh`/`height` are plain numbers, they are assumed to be in `cm`/`m`.

# Returns
A one-row `DataFrame`:
- `product::Tuple{Vararg{AbstractString}}`: product names, in table order.
- `volume::Tuple{Vararg{Float64}}`: volume assigned to each product, **in `m^3` as plain
  numbers** — a `Tuple`-typed column has no single `eltype` `removeunits`/`restoreunits`
  can strip/restore a unit from (the same reason [`StratifiedSampling`](@ref)'s `areah`
  column is already unitless), so this is stored already-stripped instead of silently
  carrying a hidden `Unitful` type through the round trip.
- `logs::Tuple{Vararg{Int}}`: number of logs cut for each product.
- `totalvolume::Vol`: volume assigned across every product (`Unitful`, `m^3`).
- `totallogs::Int`: total logs cut across every product.

# Technical description
A **deterministic greedy priority scan**, not a value/revenue optimizer: for
each product in table order, cut the longest permitted log until the stem
narrows below that product's `sed`, then move to the next product for the
remainder of the stem (see `_cutlogs`). This mirrors the published `timbeR`
algorithm exactly, generalized from 4 hand-duplicated per-model
implementations to one model-agnostic loop.

# Examples
```julia-repl
julia> using DataFrames

julia> dbh = [20.0, 20.0, 20.0, 30.0, 30.0, 30.0]; height = [18.0, 18.0, 18.0, 22.0, 22.0, 22.0];

julia> hi = [0.3, 6.0, 14.0, 0.3, 8.0, 18.0]; di = [22.4, 15.8, 6.1, 33.6, 24.2, 8.7];

julia> ft = fit(Kozak1969(), dbh, height, hi, di);

julia> products = DataFrame(
         name=["Sawlog", "Pulpwood"],
         sed=[18.0, 8.0],
         minlength=[2.5, 2.0],
         maxlength=[4.0, 3.0],
         kerf=[0.03, 0.03],
       );

julia> logassortment(ft, 30.0u"cm", 22.0u"m", products)
1×5 DataFrame
 Row │ product                     volume                   logs    totalvolume  totallogs
     │ Tuple…                      Tuple…                   Tuple…  Vol…         Int64
─────┼───────────────────────────────────────────────────────────────────────────────────
   1 │ ("Sawlog", "Pulpwood")      (0.284872, 0.0908652)     (2, 3)  0.375737 m^3          5
```
"""
function logassortment(fit::TaperFit, dbh::Len, height::Len, products::AbstractDataFrame; stumpheight::Len=0.0 * HUNIT)
  required = (:name, :sed, :minlength, :maxlength, :kerf)
  issubset(required, propertynames(products)) ||
    throw(ArgumentError("products must have columns $required."))

  name = String.(products.name)
  sed = ustrip.(DUNIT, withunit(products.sed, DUNIT))
  minlen = ustrip.(HUNIT, withunit(products.minlength, HUNIT))
  maxlen = ustrip.(HUNIT, withunit(products.maxlength, HUNIT))
  kerf = ustrip.(HUNIT, withunit(products.kerf, HUNIT))

  dbhval, heightval = ustrip(DUNIT, dbh), ustrip(HUNIT, height)
  stumpval = ustrip(HUNIT, stumpheight)
  (0 <= stumpval < heightval) || throw(DomainError(stumpheight, "stumpheight must satisfy 0 <= stumpheight < height."))

  volumes, logcounts = _cutlogs(fit, dbhval, heightval, stumpval, sed, minlen, maxlen, kerf)

  return DataFrame(
    product=Tuple(name), volume=Tuple(volumes), logs=Tuple(logcounts),
    totalvolume=sum(volumes) * VUNIT, totallogs=sum(logcounts),
  )
end

logassortment(fit::TaperFit, dbh::Real, height::Real, products::AbstractDataFrame; stumpheight::Real=0.0) =
  logassortment(fit, dbh * DUNIT, height * HUNIT, products; stumpheight=stumpheight * HUNIT)
