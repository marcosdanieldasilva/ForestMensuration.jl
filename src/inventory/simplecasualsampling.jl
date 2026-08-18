# A plot occupying a full reference area (1 ha in the metric system, 1 ac in the
# imperial one) is treated as effectively infinite regardless of the sampling fraction,
# matching the convention used by the original tool this design was ported from.
_isreferencearea(a::Area) = a == 1u"ha" || a == 1u"ac"

function sampling(::Type{SimpleCasualSampling}, volume::AbstractVector{<:Vol}, plot_area::Area, total_area::Area;
  e::Real=10, α::Real=0.95)

  n = length(volume)
  n > 1 || throw(ArgumentError("At least two sampled plots are required."))

  N = round(Int, ustrip(uconvert(NoUnits, total_area / plot_area)))
  f = 1 - n / N
  x̅ = mean(volume)
  cv = std(volume) / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N) || _isreferencearea(plot_area)
  if infinite
    population = "infinite"
    s²x̅ = var(volume) / n
    requiredplots = _requiredsamplesize(cv^2, n, e, α)
  else
    population = "finite"
    s²x̅ = (var(volume) / n) * f
    requiredplots = _requiredsamplesize(cv^2, n, e, α, N)
  end
  sx̅ = sqrt(s²x̅)
  missingplots = n > requiredplots ? 0 : requiredplots - n

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevolume = x̅ / plot_area
  totalvolume = x̅ * N
  ciupper = totalvolume + N * absoluteerror
  cilower = totalvolume - N * absoluteerror

  return DataFrame(
    vm=x̅, cv=round(cv, digits=2), s2m=s²x̅, se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vha=hectarevolume, vtotal=totalvolume,
    cilower=cilower, ciupper=ciupper, pop=population, f=round(f, digits=3),
    n=n, nreq=requiredplots, nmiss=missingplots, N=N, area=total_area,
  )
end

sampling(::Type{SimpleCasualSampling}, volume::AbstractVector{<:Real}, plot_area::Real, total_area::Real; kwargs...) =
  sampling(SimpleCasualSampling, volume * VUNIT, plot_area * AUNIT, total_area * AUNIT; kwargs...)
