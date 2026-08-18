function sampling(::Type{DoubleSampling}, volume1::AbstractVector{<:Vol}, volume2::AbstractVector{<:Union{Missing,Vol}}, plot_area::Area,
  N::Integer; α::Real=0.95)

  length(volume1) == length(volume2) || throw(DimensionMismatch("volume1 and volume2 must have the same length."))
  permanentidx = findall(!ismissing, volume2)
  m = length(permanentidx)
  n = length(volume1)
  m >= 3 || throw(ArgumentError("At least three permanent (remeasured) plots are required."))

  o1 = _occasionmean(volume1, N)
  occasion1 = _occasiontable(o1, N, plot_area, α)
  vunit = unit(o1.x̅)

  x1perm = ustrip.(volume1[permanentidx])
  y2perm = ustrip.(volume2[permanentidx])
  x̄1perm, ȳ2perm = mean(x1perm), mean(y2perm)
  x̄1all = ustrip(o1.x̅)

  b = _olsslope(y2perm, x1perm)
  sstY = sum(abs2, y2perm .- ȳ2perm)
  ssxy = sum((x1perm .- x̄1perm) .* (y2perm .- ȳ2perm))
  syx² = (sstY - b * ssxy) / (m - 2)
  sy² = sstY / (m - 1)

  yreg = ȳ2perm + b * (x̄1all - x̄1perm)
  meanvar = syx² / m + (sy² - syx²) / n
  yregq, sx̅ = yreg * vunit, sqrt(meanvar) * vunit

  t = -quantile(TDist(m - 1), (1 - α) / 2)
  absoluteerror = t * sx̅
  relativeerror = absoluteerror / yregq * 100
  hectarevolume = yregq / plot_area
  totalvolume = yregq * N

  occasion2 = DataFrame(
    vm=yregq, b=round(b, digits=4), s2yx=syx² * vunit^2, se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vha=hectarevolume, vtotal=totalvolume,
    cilower=totalvolume - N * absoluteerror, ciupper=totalvolume + N * absoluteerror,
    m=m, ntemp=n - m, N=N,
  )

  growth = yregq - o1.x̅
  change = _changetable(growth, meanvar, m - 1, N, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function sampling(::Type{DoubleSampling}, volume1::AbstractVector{<:Real}, volume2::AbstractVector{<:Union{Missing,Real}}, plot_area::Real,
  N::Integer; kwargs...)
  sampling(DoubleSampling, volume1 * VUNIT, _asvolume(volume2), plot_area * AUNIT, N; kwargs...)
end
