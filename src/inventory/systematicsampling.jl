# Successive differences between consecutive plots along a systematic line, skipping the
# pair that straddles two different lines when a `line` grouping is given -- generic to
# any number of lines (transects), reducing to a single running difference when there is
# only one.
function _successivedifferences(vol::AbstractVector{<:Vol}, line::Union{AbstractVector,Nothing})
  n = length(vol)
  keep = isnothing(line) ? trues(n - 1) : [line[i] == line[i+1] for i in 1:n-1]
  return [vol[i] - vol[i+1] for i in 1:n-1 if keep[i]]
end

function sampling(::Type{SystematicSampling}, volume::AbstractVector{<:Vol}, plot_area::Area, total_area::Area;
  line::Union{AbstractVector,Nothing}=nothing, e::Real=10, α::Real=0.95)

  n = length(volume)
  n > 1 || throw(ArgumentError("At least two sampled plots are required."))
  isnothing(line) || length(line) == n || throw(DimensionMismatch("line must have the same length as volume."))

  N = round(Int, ustrip(uconvert(NoUnits, total_area / plot_area)))
  k = isnothing(line) ? 1 : length(unique(line))
  x̅ = mean(volume)
  diffs = _successivedifferences(volume, line)

  infinite = _isinfinitepopulation(n, N)
  f = infinite ? 1.0 : 1 - n / N
  population = infinite ? "infinite" : "finite"
  s²x̅ = sum(abs2, diffs) / (2n * (n - k)) * f
  sx̅ = sqrt(s²x̅)
  cv = sx̅ / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  hectarevolume = x̅ / plot_area
  totalvolume = x̅ * N
  ciupper = totalvolume + N * absoluteerror
  cilower = totalvolume - N * absoluteerror

  return DataFrame(
    vm=x̅, cv=round(cv, digits=2), s2m=s²x̅, se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), ereq=e, vha=hectarevolume, vtotal=totalvolume,
    cilower=cilower, ciupper=ciupper, pop=population, f=round(f, digits=3),
    k=k, n=n, N=N, area=total_area,
  )
end

function sampling(::Type{SystematicSampling}, volume::AbstractVector{<:Real}, plot_area::Real, total_area::Real; kwargs...)
  sampling(SystematicSampling, volume * VUNIT, plot_area * AUNIT, total_area * AUNIT; kwargs...)
end
