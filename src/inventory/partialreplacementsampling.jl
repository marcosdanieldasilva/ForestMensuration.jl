function sampling(::Type{PartialReplacementSampling}, volume1::AbstractVector{<:Union{Missing,Vol}}, volume2::AbstractVector{<:Union{Missing,Vol}},
  plot_area::Area, N::Integer; α::Real=0.95)

  length(volume1) == length(volume2) || throw(DimensionMismatch("volume1 and volume2 must have the same length."))

  idx = eachindex(volume1)
  matched = findall(i -> !ismissing(volume1[i]) && !ismissing(volume2[i]), idx)
  unmatched = findall(i -> !ismissing(volume1[i]) && ismissing(volume2[i]), idx)
  newtemp = findall(i -> ismissing(volume1[i]) && !ismissing(volume2[i]), idx)
  m, u, v = length(matched), length(unmatched), length(newtemp)
  m >= 3 || throw(ArgumentError("At least three matched (remeasured on both occasions) plots are required."))
  v >= 2 || throw(ArgumentError("At least two new temporary (occasion-2-only) plots are required."))

  o1 = _occasionmean(collect(skipmissing(volume1)), N)
  occasion1 = _occasiontable(o1, N, plot_area, α)
  vunit = unit(o1.x̅)

  x1m, y2m = ustrip.(volume1[matched]), ustrip.(volume2[matched])
  x1u = ustrip.(volume1[unmatched])
  y2v = ustrip.(volume2[newtemp])
  n1 = u + m
  x̅u = isempty(x1u) ? zero(x1m[1]) : mean(x1u)
  x̅m, ȳm = mean(x1m), mean(y2m)
  x̅combined = (u * x̅u + m * x̅m) / n1
  Sym² = var(y2m)

  b = _olsslope(y2m, x1m)
  sstY = sum(abs2, y2m .- ȳm)
  ssxy = sum((x1m .- x̅m) .* (y2m .- ȳm))
  syx² = (sstY - b * ssxy) / (m - 2)
  yreg = ȳm + b * (x̅combined - x̅m)
  varyreg = syx² / m + (Sym² - syx²) / n1

  ȳv = mean(y2v)
  varyv = var(y2v) / v * (1 - v / N)

  c = varyv / (varyreg + varyv)
  yspr = c * yreg + (1 - c) * ȳv
  varyspr = c^2 * varyreg + (1 - c)^2 * varyv

  ysprq, sx̅ = yspr * vunit, sqrt(varyspr) * vunit
  t = -quantile(TDist(m + v - 2), (1 - α) / 2)
  absoluteerror = t * sx̅
  relativeerror = absoluteerror / ysprq * 100
  hectarevolume = ysprq / plot_area
  totalvolume = ysprq * N

  occasion2 = DataFrame(
    vm=ysprq, vreg=yreg * vunit, vtemp=ȳv * vunit, c=round(c, digits=4), b=round(b, digits=4),
    se=sx̅, abserr=absoluteerror, relerr=round(relativeerror, digits=2), vha=hectarevolume,
    vtotal=totalvolume, cilower=totalvolume - N * absoluteerror, ciupper=totalvolume + N * absoluteerror,
    m=m, u=u, v=v, N=N,
  )

  growth = ysprq - x̅combined * vunit
  change = _changetable(growth, varyspr, m + v - 2, N, α)

  return SamplingReport((; occasion1, occasion2, change))
end

function sampling(::Type{PartialReplacementSampling}, volume1::AbstractVector{<:Union{Missing,Real}}, volume2::AbstractVector{<:Union{Missing,Real}},
  plot_area::Real, N::Integer; kwargs...)
  sampling(PartialReplacementSampling, _asvolume(volume1), _asvolume(volume2), plot_area * AUNIT, N; kwargs...)
end
