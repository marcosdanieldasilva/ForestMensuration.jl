# Per-point sample: sums the individual tree contributions -- basal area/ha and
# volume/ha -- for every point that had at least one "in" tree. Points with zero "in"
# trees never appear as rows in `data` and are added back as explicit zero observations
# by the caller (see `nzero` below) -- omitting them would systematically inflate the
# mean, a well-known pitfall when processing angle-count data.
function _pointtable(point::Symbol, ef::AbstractVector{<:Quantity}, g::AbstractVector{<:Area},
  vol::AbstractVector{<:Vol}, data::AbstractDataFrame)

  working = DataFrame(point=data[!, point], ef=ef, g=g, vol=vol)
  return combine(groupby(working, :point, sort=true)) do df
    (ntrees=nrow(df), Gha=sum(df.ef .* df.g), vha=sum(df.ef .* df.vol))
  end
end

function sampling(::Type{HorizontalPointSampling}, point::Symbol, diameter::Symbol, volume::Symbol, baf::Real,
  npoints::Integer, point_area::Area, total_area::Area, data::AbstractDataFrame; e::Real=10, α::Real=0.95)

  baf > 0 || throw(ArgumentError("baf must be a positive value."))

  d = _asdiameter(data[!, diameter])
  vol = _asvolume(data[!, volume])
  g = basalarea.(d)
  bafunit = unit(g[1]) / referencearea(unit(d[1]))
  bafQ = baf * bafunit
  ef = bafQ ./ g

  pointtable = _pointtable(point, ef, g, vol, data)
  ncounted = nrow(pointtable)
  nzero = npoints - ncounted
  nzero >= 0 || throw(ArgumentError(
    "npoints ($npoints) is smaller than the number of distinct points with at least " *
    "one tree in data[!, point] ($ncounted)."))

  n = npoints
  n > 1 || throw(ArgumentError("At least two sampled points are required."))

  vhavalues = nzero == 0 ? pointtable.vha : vcat(pointtable.vha, fill(0.0, nzero) .* unit(pointtable.vha[1]))
  Ghavalues = nzero == 0 ? pointtable.Gha : vcat(pointtable.Gha, fill(0.0, nzero) .* unit(pointtable.Gha[1]))

  N = round(Int, ustrip(uconvert(NoUnits, total_area / point_area)))
  f = 1 - n / N
  x̅ = mean(vhavalues)
  cv = std(vhavalues) / x̅ * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)

  infinite = _isinfinitepopulation(n, N)
  if infinite
    population = "infinite"
    s²x̅ = var(vhavalues) / n
    requiredpoints = _requiredsamplesize(cv^2, n, e, α)
  else
    population = "finite"
    s²x̅ = (var(vhavalues) / n) * f
    requiredpoints = _requiredsamplesize(cv^2, n, e, α, N)
  end
  sx̅ = sqrt(s²x̅)
  missingpoints = n > requiredpoints ? 0 : requiredpoints - n

  absoluteerror = Ttab * sx̅
  relativeerror = absoluteerror / x̅ * 100
  totalvolume = x̅ * total_area
  ciupper = totalvolume + total_area * absoluteerror
  cilower = totalvolume - total_area * absoluteerror
  Gha = mean(Ghavalues)
  Gtotal = Gha * total_area

  resulttable = DataFrame(
    vha=x̅, cv=round(cv, digits=2), s2m=s²x̅, se=sx̅, abserr=absoluteerror,
    relerr=round(relativeerror, digits=2), vtotal=totalvolume, cilower=cilower, ciupper=ciupper,
    Gha=Gha, Gtotal=Gtotal, pop=population, f=round(f, digits=3),
    n=n, nreq=requiredpoints, nmiss=missingpoints, N=N, area=total_area,
  )

  return SamplingReport((; pointTable=pointtable, resultTable=resulttable))
end

function sampling(::Type{HorizontalPointSampling}, point::Symbol, diameter::Symbol, volume::Symbol, baf::Real,
  npoints::Integer, point_area::Real, total_area::Real, data::AbstractDataFrame; kwargs...)
  sampling(HorizontalPointSampling, point, diameter, volume, baf, npoints, point_area * AUNIT, total_area * AUNIT, data; kwargs...)
end
