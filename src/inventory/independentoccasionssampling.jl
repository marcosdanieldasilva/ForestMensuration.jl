function sampling(::Type{IndependentOccasionsSampling}, volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Vol}, plot_area::Area,
  N1::Integer, N2::Integer; α::Real=0.95)

  o1 = _occasionmean(volume1, N1)
  o2 = _occasionmean(volume2, N2)
  occasion1 = _occasiontable(o1, N1, plot_area, α)
  occasion2 = _occasiontable(o2, N2, plot_area, α)

  growth = o2.x̅ - o1.x̅
  growthvar = ustrip(o1.meanvar) + ustrip(o2.meanvar)
  df = (o1.n - 1) + (o2.n - 1) - 1
  change = _changetable(growth, growthvar, df, N2, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function sampling(::Type{IndependentOccasionsSampling}, volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Real}, plot_area::Real,
  N1::Integer, N2::Integer; kwargs...)
  sampling(IndependentOccasionsSampling, volume1 * VUNIT, volume2 * VUNIT, plot_area * AUNIT, N1, N2; kwargs...)
end
