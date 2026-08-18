# Internal helpers shared by every sampling design in this submodule. None of these are
# exported; each public function calls into them with the specific quantities its own
# variance formula produces.

# Normalizes a DataFrame column to a `Vol` vector regardless of whether it already
# carries units: used by every grouped (Symbol + DataFrame) sampling design, since a
# column's element type isn't known until the DataFrame is indexed at runtime, unlike a
# directly-dispatched `AbstractVector` argument.
_asvolume(v::AbstractVector{<:Vol}) = v
_asvolume(v::AbstractVector{<:Real}) = v .* VUNIT
_asvolume(v::AbstractVector{<:Union{Missing,Vol}}) = v
_asvolume(v::AbstractVector{<:Union{Missing,Real}}) = [ismissing(x) ? missing : x * VUNIT for x in v]

# Same normalization as `_asvolume`, for diameter columns -- used by
# `HorizontalPointSampling`, the only design in this submodule that works from
# individual-tree data rather than pre-aggregated plot volumes.
_asdiameter(d::AbstractVector{<:Len}) = d
_asdiameter(d::AbstractVector{<:Real}) = d .* DUNIT
_asdiameter(d::AbstractVector{<:Union{Missing,Len}}) = d
_asdiameter(d::AbstractVector{<:Union{Missing,Real}}) = [ismissing(x) ? missing : x * DUNIT for x in d]

# Required sample size for an infinite population, given a t-value, a *squared*
# variance-like term (`varterm`, e.g. `cv^2` or a variance in volume² units -- not `cv` or
# a standard deviation) and an admissible error (`e`) on the matching linear scale (both
# relative, e.g. percent of the mean, or both absolute, e.g. volume units).
_infinitesamplesize(t::Real, varterm::Real, e::Real) = (t^2 * varterm) / e^2

# Required sample size for a finite population of size N, on the same scale as above.
# `correctionvarterm` lets the finite-population correction term use a different
# variance-like quantity than the numerator -- needed by `TwoStageSampling`, where the
# correction reflects the population size of secondary units while the numerator reflects
# how many were actually measured; it defaults to `varterm` for every other design, where
# the two coincide.
_finitesamplesize(t::Real, varterm::Real, e::Real, N::Real, correctionvarterm::Real=varterm) =
  (t^2 * varterm) / (e^2 + (t^2 * correctionvarterm) / N)

# Iteratively solves for the required sample size. The t quantile depends on the degrees
# of freedom of the very sample size it is used to compute, so the pair is refined
# together until the estimate stabilizes (matching the fixed-point iteration used
# throughout the legacy ForestInventory.jl sampling routines). `N=nothing` treats the
# population as infinite; passing a population size switches to the finite-population
# formula. Shared by every sampling design below instead of each reimplementing the loop.
function _requiredsamplesize(varterm::Real, n::Real, e::Real, α::Real, N::Union{Real,Nothing}=nothing;
  correctionvarterm::Real=varterm, max_iterations::Int=1000)
  samplesize(t::Real) = isnothing(N) ? _infinitesamplesize(t, varterm, e) : _finitesamplesize(t, varterm, e, N, correctionvarterm)
  for _ in 1:max_iterations
    t1 = -quantile(TDist(n - 1), (1 - α) / 2)
    iter1 = samplesize(t1)
    t2 = -quantile(TDist(iter1 - 1), (1 - α) / 2)
    iter2 = samplesize(t2)
    if round(n, digits=3) == round(iter1, digits=3) || round(n, digits=3) == round(iter2, digits=3)
      return ceil(Int, n)
    end
    n = iter1
  end
  error("Maximum number of iterations ($max_iterations) reached without convergence.")
end

# Whether a sample of size `n` drawn without replacement from a population of size `N`
# can be treated as effectively infinite: the standard forest-mensuration rule of thumb
# is that the finite-population correction is negligible once the sampling fraction stays
# under 2%.
_isinfinitepopulation(n::Real, N::Real) = (1 - n / N) >= 0.98

# Between/within primary-unit variance decomposition shared by every design built on
# equal-size primary sampling units (clusters, systematic starts, or first-stage units of
# a two-stage design): `means`/`vars` are the per-primary-unit mean and internal variance
# (one entry per sampled primary unit), `M` is the number of secondary units per primary
# unit (assumed equal), and `N` is the number of possible primary units in the population.
#
# Returns a named tuple with the pooled within-unit variance, the between-unit variance
# (per secondary unit), their sum (the total per-secondary-unit variance), the intraclass
# correlation coefficient, and the variance of the overall mean.
function _clustervariance(means::AbstractVector{<:Real}, vars::AbstractVector{<:Real}, M::Real, N::Real)
  n = length(means)
  ȳ = mean(means)
  # a cluster of exactly one secondary unit has no internal spread to estimate (`var` of a
  # single observation is NaN, not 0); with M=1 there simply is no within-unit variance,
  # and this reduces the whole decomposition to simple random sampling of the means.
  withinvar = M == 1 ? zero(eltype(vars)) : mean(vars)  # pooled within-unit variance: Σ(M-1)sⱼ² / (n(M-1)) simplifies to mean(sⱼ²) for equal M
  msb = M * sum(abs2, means .- ȳ) / (n - 1)  # between-unit mean square
  betweenvar = (msb - withinvar) / M
  totalvar = betweenvar + withinvar
  icc = totalvar > 0 ? betweenvar / totalvar : zero(betweenvar)
  f = _isinfinitepopulation(n, N) ? one(N) : (N - n) / N
  # `f` scales the whole between+within bracket, not just the between term: with a full
  # census inside every sampled cluster there is a single source of sampling error (which
  # clusters were drawn), and `f` is that one correction applied to the combined estimate.
  meanvar = f * (betweenvar / n + withinvar / (n * M))
  varterm = betweenvar + withinvar / M  # feeds `_requiredsamplesize` directly
  return (; withinvar, betweenvar, totalvar, icc, meanvar, varterm)
end

# Simple-random-sampling point estimate and variance of the mean for one occasion of a
# successive-occasions design (independent/complete-replacement/partial-replacement/
# double sampling), reused instead of recomputing it for every occasion in every design.
function _occasionmean(vol::AbstractVector{<:Vol}, N::Real)
  n = length(vol)
  x̅ = mean(vol)
  s² = var(vol)
  meanvar = (s² / n) * (1 - n / N)
  return (; n, x̅, s², meanvar)
end

# Builds one occasion's result table (mean, standard error, confidence interval, total
# volume) as a single-row wide `DataFrame`, shared by every successive-occasions design
# below. Column names are short, ASCII, and unit-homogeneous per column -- exactly the
# shape `removeunits`/`restoreunits` (ForestFoundations.jl) expect, and the same
# convention as `dmetrics`/`standmetrics` elsewhere in this package.
function _occasiontable(stats::NamedTuple, N::Real, plot_area::Area, α::Real)
  t = -quantile(TDist(stats.n - 1), (1 - α) / 2)
  se = sqrt(stats.meanvar)
  abserr = t * se
  relerr = abserr / stats.x̅ * 100
  vtotal = stats.x̅ * N
  vha = stats.x̅ / plot_area
  return DataFrame(
    vm=stats.x̅, se=se, abserr=abserr, relerr=round(relerr, digits=2), vha=vha,
    vtotal=vtotal, cilower=vtotal - N * abserr, ciupper=vtotal + N * abserr,
    n=stats.n, N=N,
  )
end

# Builds the change-between-occasions result table (mean growth, standard error,
# confidence interval, total growth) as a single-row wide `DataFrame`, shared by every
# successive-occasions design below.
function _changetable(growth::Vol, growthvar::Real, df::Real, N::Real, α::Real)
  t = -quantile(TDist(df), (1 - α) / 2)
  se = sqrt(growthvar) * unit(growth)
  abserr = t * se
  relerr = abserr / growth * 100
  gtotal = growth * N
  return DataFrame(
    gm=growth, se=se, abserr=abserr, relerr=round(relerr, digits=2),
    gtotal=gtotal, cilower=gtotal - N * abserr, ciupper=gtotal + N * abserr,
  )
end

# Ordinary least squares slope of `y` on `x` via the normal equations, solved through the
# same transposed-design-matrix approach as the ANOVA helper in `stratifiedsampling.jl`
# (`X'X`, `X'y`, Cholesky solve) rather than a dedicated regression model — this module
# only ever needs the slope of a single predictor, used as the regression coefficient in
# the double-sampling and partial-replacement estimators.
function _olsslope(y::AbstractVector{<:Real}, x::AbstractVector{<:Real})
  length(y) == length(x) || throw(DimensionMismatch("y and x must have the same length."))
  X = [ones(length(x)) x]
  β = X'y
  chol = cholesky!(Symmetric(X'X))
  ldiv!(chol, β)
  return β[2]
end
