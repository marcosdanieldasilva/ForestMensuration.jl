# Two-stage sampling's variance decomposition, unlike one-stage cluster sampling's
# `_clustervariance`, has two independent sources of sampling error: which `n` (of `N`)
# primary units were drawn, and which `m` (of `M`) secondary units were measured within
# each of them. Cochran (1977), eq. 10.4, equal-size primaries, SRS at both stages.
function _twostagevariance(means::AbstractVector{<:Real}, vars::AbstractVector{<:Real}, m::Real, M::Real, N::Real)
  n = length(means)
  ȳ = mean(means)
  withinvar = mean(vars)
  msb = m * sum(abs2, means .- ȳ) / (n - 1)
  betweenvar = (msb - withinvar) / m
  f1, f2 = n / N, m / M
  meanvar = (1 - f1) / n * betweenvar + f1 * (1 - f2) / (n * m) * withinvar
  varterm = betweenvar + withinvar / m
  correctionvarterm = betweenvar + withinvar / M
  return (; withinvar, betweenvar, meanvar, varterm, correctionvarterm)
end

function sampling(::Type{TwoStageSampling}, primary::Symbol, volume::Symbol, plot_area::Area, N::Integer, M::Integer,
  data::AbstractDataFrame; e::Real=10, α::Real=0.95)

  vol = _asvolume(data[!, volume])
  table = _clustertable(primary, vol, data)
  m = _equalclustersize(table.n)
  n = nrow(table)
  m <= M || throw(ArgumentError("The number of measured secondary units per primary unit (m=$m) cannot exceed the population size M=$M."))

  x̅ = mean(table.x̅)
  vunit = unit(x̅)
  v2unit = vunit^2

  variance = _twostagevariance(ustrip.(table.x̅), ustrip.(table.s²), m, M, N)
  cv = sqrt(variance.withinvar + variance.betweenvar) / ustrip(x̅) * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  E = (e / 100) * x̅
  requiredprimaries = _requiredsamplesize(variance.varterm, n, ustrip(E), α, infinite ? nothing : N;
    correctionvarterm=variance.correctionvarterm)
  missingprimaries = n > requiredprimaries ? 0 : requiredprimaries - n

  population = infinite ? "infinite" : "finite"
  sx̅ = sqrt(variance.meanvar) * vunit
  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevolume = x̅ / plot_area
  totalvolume = x̅ * N * M
  ciupper = totalvolume + N * M * absoluteerror
  cilower = totalvolume - N * M * absoluteerror

  resulttable = DataFrame(
    vm=x̅, cv=round(cv, digits=2), s2w=variance.withinvar * v2unit, s2b=variance.betweenvar * v2unit,
    se=sx̅, abserr=absoluteerror, relerr=round(relativeerror, digits=2), vha=hectarevolume,
    vtotal=totalvolume, cilower=cilower, ciupper=ciupper, pop=population, f=round(1 - n / N, digits=3),
    m=m, M=M, n=n, nreq=requiredprimaries, nmiss=missingprimaries, N=N,
  )

  return SamplingReport((; primaryTable=table, resultTable=resulttable))
end

function sampling(::Type{TwoStageSampling}, primary::Symbol, volume::Symbol, plot_area::Real, N::Integer, M::Integer,
  data::AbstractDataFrame; kwargs...)
  sampling(TwoStageSampling, primary, volume, plot_area * AUNIT, N, M, data; kwargs...)
end
