# Per-stratum descriptive statistics: sample size, mean, variance, standard deviation,
# and the stratum's relative weight `p` (its share of the total population area) together
# with `p·s` and `p·s²`, the building blocks of every stratified-sampling formula below.
# Generic to any number of strata via `groupby`.
function _auxiliarytable(stratum::Symbol, vol::AbstractVector{<:Vol}, strata_area::AbstractVector{<:Area}, data::AbstractDataFrame)
  data = copy(data)
  data[!, :__volume__] = vol
  table = combine(groupby(data, stratum, sort=true)) do df
    (n=nrow(df), x̅=mean(df.__volume__), s²=var(df.__volume__), s=std(df.__volume__))
  end
  p = strata_area ./ sum(strata_area)
  table.p = p
  table.ps = p .* table.s
  table.ps² = p .* table.s²
  return table
end

# One-way ANOVA testing whether the strata means differ, via the standard sum-of-squares
# decomposition -- mathematically identical to fitting `volume ~ 1 + stratum` and reading
# off its F-test, without needing a modeling/design-matrix dependency for what is, here,
# a single categorical predictor.
function _anova(vol::AbstractVector{<:Vol}, table::AbstractDataFrame)
  n, k = length(vol), nrow(table)
  ȳ = mean(vol)
  sst = sum(abs2, vol .- ȳ)
  ssb = sum(table.n .* abs2.(table.x̅ .- ȳ))
  ssw = sst - ssb
  dfb, dfw = k - 1, n - k
  msb, msw = ssb / dfb, ssw / dfw
  # a single stratum has no between-strata comparison to make: the F-test is undefined
  # (0 numerator degrees of freedom), so the between-strata row is left blank instead.
  fstat = dfb > 0 ? msb / msw : missing
  pval = dfb > 0 ? ccdf(FDist(dfb, dfw), fstat) : missing
  DataFrame(
    "Source" => ["Between strata", "Within strata", "Total"],
    "DOF" => [dfb, dfw, n - 1],
    "SS" => [ssb, ssw, sst],
    "MS" => [dfb > 0 ? msb : missing, msw, missing],
    "F" => [fstat, missing, missing],
    "p(>F)" => [pval, missing, missing],
  )
end

# Elementwise max(required - measured, 0): how many additional plots each stratum still
# needs, never negative when a stratum already exceeds its required allocation.
_missingplots(measured::AbstractVector{<:Integer}, required::AbstractVector{<:Integer}) = max.(required .- measured, 0)

function sampling(::Type{StratifiedSampling}, stratum::Symbol, volume::Symbol, plot_area::Area, strata_area::AbstractVector{<:Area},
  data::AbstractDataFrame; e::Real=10, α::Real=0.95)

  length(strata_area) == length(unique(data[!, stratum])) ||
    throw(ArgumentError("strata_area must have one entry per unique stratum, in sorted order."))

  vol = _asvolume(data[!, volume])
  table = _auxiliarytable(stratum, vol, strata_area, data)
  anova = _anova(vol, table)

  N = round(Int, ustrip(uconvert(NoUnits, sum(strata_area) / plot_area)))
  x̅st = sum(table.p .* table.x̅)
  s²x̅st = sum(table.ps .^ 2 ./ table.n) - sum(table.ps²) / N
  sx̅st = sqrt(s²x̅st)
  cv = sx̅st / x̅st * 100
  E = (e / 100) * x̅st
  measured = table.n
  n = sum(measured)

  infinite = _isinfinitepopulation(n, N)
  requiredtotal = _requiredsamplesize(ustrip(sum(table.ps²)), n, ustrip(E), α, infinite ? nothing : N)
  required = round.(Int, requiredtotal .* table.ps ./ sum(table.ps))
  missingplots = _missingplots(measured, required)

  Nh = ustrip.(uconvert.(NoUnits, strata_area ./ plot_area))
  gh = Nh .* (Nh .- measured) ./ measured
  dfeff = sum(gh .* table.s²)^2 / sum((gh .^ 2 .* table.s² .^ 2) ./ (measured .- 1))
  Ttab = -quantile(TDist(dfeff), (1 - α) / 2)

  population = infinite ? "infinite" : "finite"
  f = 1 - n / N
  absoluteerror = Ttab * sx̅st
  relativeerror = absoluteerror / x̅st * 100
  totalvolume = x̅st * N
  ciupper = totalvolume + N * absoluteerror
  cilower = totalvolume - N * absoluteerror

  resulttable = DataFrame(
    vm=x̅st, cv=round(cv, digits=2), s2m=s²x̅st, se=sx̅st, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vtotal=totalvolume, cilower=cilower, ciupper=ciupper,
    pop=population, f=round(f, digits=3), nh=Tuple(measured), nreqh=Tuple(required),
    nmissh=Tuple(missingplots), n=n, nreq=requiredtotal, N=N, areah=Tuple(ustrip.(strata_area)),
  )

  return SamplingReport((; anova, auxiliaryTable=table, resultTable=resulttable))
end

function sampling(::Type{StratifiedSampling}, stratum::Symbol, volume::Symbol, plot_area::Real, strata_area::AbstractVector{<:Real},
  data::AbstractDataFrame; kwargs...)
  sampling(StratifiedSampling, stratum, volume, plot_area * AUNIT, strata_area * AUNIT, data; kwargs...)
end
