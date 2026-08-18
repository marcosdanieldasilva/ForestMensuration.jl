function sampling(::Type{MultistartSystematicSampling}, start::Symbol, volume::Symbol, plot_area::Area, total_area::Area,
  data::AbstractDataFrame; e::Real=10, α::Real=0.95)

  vol = _asvolume(data[!, volume])
  table = _clustertable(start, vol, data)
  M = _equalclustersize(table.n)
  n = nrow(table)

  N = round(Int, ustrip(uconvert(NoUnits, total_area / (plot_area * M))))
  x̅ = mean(table.x̅)
  vunit = unit(x̅)
  v2unit = vunit^2

  variance = _clustervariance(ustrip.(table.x̅), ustrip.(table.s²), M, N)
  cv = sqrt(variance.totalvar) / ustrip(x̅) * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  E = (e / 100) * x̅
  requiredplots = _requiredsamplesize(variance.varterm, n, ustrip(E), α, infinite ? nothing : N)
  missingplots = n > requiredplots ? 0 : requiredplots - n

  population = infinite ? "infinite" : "finite"
  sx̅ = sqrt(variance.meanvar) * vunit
  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevolume = x̅ / plot_area
  totalvolume = x̅ * M * N
  ciupper = totalvolume + M * N * absoluteerror
  cilower = totalvolume - M * N * absoluteerror

  resulttable = DataFrame(
    vm=x̅, cv=round(cv, digits=2), s2w=variance.withinvar * v2unit, s2b=variance.betweenvar * v2unit,
    s2=variance.totalvar * v2unit, icc=round(variance.icc, digits=4), se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vha=hectarevolume, vtotal=totalvolume, cilower=cilower,
    ciupper=ciupper, pop=population, f=round(1 - n / N, digits=3), M=M, n=n, nreq=requiredplots,
    nmiss=missingplots, N=N, area=total_area,
  )

  return SamplingReport((; startTable=table, resultTable=resulttable))
end

function sampling(::Type{MultistartSystematicSampling}, start::Symbol, volume::Symbol, plot_area::Real, total_area::Real,
  data::AbstractDataFrame; kwargs...)
  sampling(MultistartSystematicSampling, start, volume, plot_area * AUNIT, total_area * AUNIT, data; kwargs...)
end
