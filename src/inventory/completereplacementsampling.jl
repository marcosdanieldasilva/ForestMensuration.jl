function sampling(::Type{CompleteReplacementSampling}, volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Vol}, plot_area::Area,
  N1::Integer, N2::Integer; α::Real=0.95)

  length(volume1) == length(volume2) ||
    throw(DimensionMismatch("volume1 and volume2 must have the same length: the same plots are remeasured on both occasions."))

  o1 = _occasionmean(volume1, N1)
  o2 = _occasionmean(volume2, N2)
  occasion1 = _occasiontable(o1, N1, plot_area, α)
  occasion2 = _occasiontable(o2, N2, plot_area, α)

  n = o1.n
  growth = o2.x̅ - o1.x̅
  growthcovariance = ustrip(cov(volume1, volume2))
  # with a small sample and strong occasion-to-occasion correlation, the covariance
  # estimate can occasionally overshoot enough to make this formula dip below zero;
  # clamped at zero, since a variance estimate cannot legitimately be negative.
  growthvar = max(ustrip(o1.meanvar) + ustrip(o2.meanvar) - 2 * growthcovariance / n, 0.0)
  df = 2 * (n - 1) - 1
  change = _changetable(growth, growthvar, df, N2, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function sampling(::Type{CompleteReplacementSampling}, volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Real}, plot_area::Real,
  N1::Integer, N2::Integer; kwargs...)
  sampling(CompleteReplacementSampling, volume1 * VUNIT, volume2 * VUNIT, plot_area * AUNIT, N1, N2; kwargs...)
end
