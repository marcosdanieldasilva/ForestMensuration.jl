# ==============================================================================
# SITE CLASSIFICATION (Guide-Curve / Delta method)
# ==============================================================================
# Anamorphic guide-curve method: shift each plot's trajectory to the index
# age by the model-implied difference in height between the current age and
# the index age. Works for any AllometricModel fitted with age as the sole
# continuous regressor, linear or transformed on either side — it does NOT
# extend to genuinely nonlinear (non-linear-in-parameters) growth models;
# that needs the Algebraic Difference Approach, out of scope here.
#
# Unit policy: `siteClassification`/`hdomClassification` return values sit
# directly on the response's (height) scale, so — like `predict` — they carry
# the response's unit when the model was fit on unitful data.
# `dataAge`/`site` inputs are matched against the model's stored units the
# same way `predict`'s `data` argument is. `siteTable` builds a wide summary
# table rather than a per-observation response value, so — like
# `criteriaTable`/`metrics` in ForestModeling.jl — it stays on the plain
# numeric (fit-unit) scale throughout, never `Unitful`-tagged.
#
# `matchunits`/`_restoreunit`/`_matchcol`/`_stripcol` are ForestModeling.jl
# internals (not exported — that package's public surface stays pure
# statistics, with no forestry-specific site/age semantics of its own),
# qualified here rather than re-exported just to serve this one caller.
# Automatic class breadth/centers reuse `_sturges`/`_amplitude`/
# `_classbreadth`/`_classcenter` from `distributiontables.jl` directly — the
# same "nice round number" convention `frequencytable`/`diametrictable`
# already use, not a second, differently-rounded reimplementation.

"""
    calculateΔ(model::AllometricModel, dataAge::AbstractDataFrame, indexAge::Real)

Model-implied difference in the (transformed-scale) response between the
index age and each observation's current age — the core of the delta method.
`dataAge` is assumed already unit-matched to `model` (see
`ForestModeling.matchunits`) — internal helper, not called directly with raw
user data.
"""
function calculateΔ(model::AllometricModel, dataAge::AbstractDataFrame, indexAge::Real)
  (yName, xName, qNames...) = propertynames(model.data)
  dataIndexAge = deepcopy(dataAge[!, [yName, xName, qNames...]])
  dataIndexAge[!, xName] .= indexAge

  mmAge = modelmatrix(model.formula.rhs.terms[2:end], dataAge)
  mmIndexAge = modelmatrix(model.formula.rhs.terms[2:end], dataIndexAge)
  β = model.β[2:end]'

  Δ = mmIndexAge .- mmAge
  return sum(β .* Δ, dims=2)[:]
end

"""
    siteClassification(model::AllometricModel, dataAge::AbstractDataFrame, indexAge::Real)
    siteClassification(model::AllometricModel, indexAge::Real)

Site index for each observation: the dominant height each plot is expected
to reach at `indexAge`, given its currently observed height/age and the
fitted growth model. The second method reuses the data the model was fit on.
Carries the response's `Unitful` unit when the model was fit on unitful data.

# Examples
```julia-repl
julia> using DataFrames

julia> ageData = DataFrame(idade=[3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0],
                            hdom=[10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2, 22.0, 23.1, 24.0, 25.0]u"m");

julia> msite = fit(AllometricModel, @formula(log(hdom) ~ 1 + idade^-1), ageData);

julia> siteClassification(msite, 10.0)   # site index at index age 10, carries hdom's unit (m)
14-element Vector{Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}:
 22.8 m
  ⋮
 22.1 m
```
"""
function siteClassification(model::AllometricModel, dataAge::AbstractDataFrame, indexAge::Real)
  indexAge <= 0 && throw(DomainError(indexAge, "index age must be positive"))
  dataAge = DataFrame(ForestModeling.matchunits(dataAge, model.units))
  Δ = calculateΔ(model, dataAge, indexAge)
  y = model.formula.lhs
  site = modelcols(y, dataAge) .+ Δ   # T(site) — still on the transformed (fitting) scale
  if isa(y, FunctionTerm)
    x = length(y.args) > 1 ? dataAge[!, y.args[2].sym] : nothing
    predictBiasCorrected!(site, x, y, model.σ²)   # back-transform in place to the original scale
  end
  return ForestModeling._restoreunit(round.(site, digits=1), model.units[1])
end
siteClassification(model::AllometricModel, indexAge::Real) = siteClassification(model, DataFrame(model.data), indexAge)

"""
    hdomClassification(model::AllometricModel, dataAge::AbstractDataFrame, indexAge::Real, site::AbstractVector)

Inverse of [`siteClassification`](@ref): given a site index and ages,
predicts dominant height — used to forecast stand development per site
class. `site` may be a plain-number vector or a `Unitful.Quantity` vector
compatible with the model's response unit (e.g. the direct output of
[`siteClassification`](@ref) — `hdomClassification(model, dados, idade,
siteClassification(model, dados, idade))` is the round-trip sanity check
this signature is meant to support).

# Examples
```julia-repl
julia> hdomClassification(msite, ageData, 10.0, siteClassification(msite, ageData, 10.0))   # ≈ ageData.hdom
14-element Vector{Quantity{Float64, 𝐋, Unitful.FreeUnits{(m,), 𝐋, nothing}}}:
 10.3 m
  ⋮
 25.2 m
```
"""
function hdomClassification(model::AllometricModel, dataAge::AbstractDataFrame, indexAge::Real, site::AbstractVector)
  indexAge <= 0 && throw(DomainError(indexAge, "index age must be positive"))
  siteValues = ForestModeling._matchcol(site, model.units[1])
  any(<(0), siteValues) && throw(DomainError(siteValues, "site values must be positive"))
  dataAge = DataFrame(ForestModeling.matchunits(dataAge, model.units))
  Δ = calculateΔ(model, dataAge, indexAge)
  (yName, xName, qNames...) = propertynames(model.data)
  siteData = deepcopy(dataAge[!, [yName, xName, qNames...]])
  siteData[!, yName] .= siteValues
  siteData[!, xName] .= indexAge
  y = model.formula.lhs
  hdom = modelcols(y, siteData) .- Δ   # T(hdom at dataAge's own ages)
  if isa(y, FunctionTerm)
    x = length(y.args) > 1 ? dataAge[!, y.args[2].sym] : nothing
    predictBiasCorrected!(hdom, x, y, model.σ²)   # back-transform in place to the original scale
  end
  return ForestModeling._restoreunit(round.(hdom, digits=1), model.units[1])
end

"""
    siteTable(model::AllometricModel, indexAge::Real, hi::Real)
    siteTable(model::AllometricModel, indexAge::Real)

Table of predicted dominant heights per age and site class (automatic
class breadth via [`frequencytable`](@ref)'s own Sturges'-rule-based
"nice number" rounding when `hi` is not given). Plotting is intentionally
left out — this returns the `DataFrame` only. Like `criteriaTable`/`metrics`
in ForestModeling.jl, this is a summary table rather than a per-observation
response value, so it stays on the plain numeric (fit-unit) scale even when
`model` was fit on unitful data.

# Examples

`msite` here is the same model fit on unitful `hdom` (meters) from
[`siteClassification`](@ref)'s example — the table's values stay plain `Float64` in the
fitting unit (meters) regardless, per the unit policy above:

```julia-repl
julia> siteTable(msite, 10.0)   # class breadth picked automatically
14×8 DataFrame
 Row │ idade    S_17.5   S_18.5   S_19.5   S_20.5   S_21.5   S_22.5   S_23.5
─────┼─────────────────────────────────────────────────────────────────────
   1 │     3.0      7.9      8.4      8.8      9.3      9.7     10.2     10.6
  ⋮  │    ⋮        ⋮        ⋮        ⋮        ⋮        ⋮        ⋮        ⋮

julia> siteTable(msite, 10.0, 2.0)   # fixed class breadth of 2.0, same convention as frequencytable(x, hi)
14×4 DataFrame
 Row │ idade    S_19.0   S_21.0   S_23.0
─────┼───────────────────────────────────
   1 │     3.0      8.6      9.5     10.4
  ⋮  │    ⋮        ⋮        ⋮        ⋮
```
"""
function siteTable(model::AllometricModel, indexAge::Real, hi::Real)
  hi <= 0 && throw(DomainError(hi, "class breadth must be positive"))
  (hd, age, q...) = propertynames(model.data)
  site = ForestModeling._stripcol(siteClassification(model, indexAge))
  sites = sort(unique(_classcenter.(site, hi)))
  ages = sort(unique(model.data[age]))

  repeatedAges = repeat(ages, outer=length(sites))
  repeatedSites = vcat([fill(s, length(ages)) for s in sites]...)
  scaffold = DataFrame(age => repeatedAges, hd => repeatedSites)
  hdomPredicted = ForestModeling._stripcol(hdomClassification(model, scaffold, indexAge, repeatedSites))
  insertcols!(scaffold, :site => hdomPredicted, makeunique=true)

  table = unstack(scaffold, propertynames(scaffold)...) |> dropmissing!
  rename!(table, [age; [Symbol("S_$s") for s in names(table)[2:end]]])
  return table
end
function siteTable(model::AllometricModel, indexAge::Real)
  site = ForestModeling._stripcol(siteClassification(model, indexAge))
  hi = _classbreadth(_amplitude(site), _sturges(length(site)))
  return siteTable(model, indexAge, hi)
end
