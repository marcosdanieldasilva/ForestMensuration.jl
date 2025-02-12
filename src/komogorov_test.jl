export ks_one_sample_test, ks_adjusted_one_sample_test, ks_two_sample_test, p_value

abstract type KolmogorovSmirnovTest end

struct KSOneSampleTest <: KolmogorovSmirnovTest
  α::Float64
  n_elements::Int
  Dplus::Float64
  Dminus::Float64
  Dn::Float64
  Dcrit::Float64
  distribution::UnivariateDistribution
  test_type::Symbol
end

name_of_test(::KSOneSampleTest) = "Kolmogorov-Smirnov One-Sample Test"

function show_params(io::IO, x::KSOneSampleTest)
  println(io, "  Number of observations:  $(x.n_elements)")
  println(io, "  Test parameters:         $(x.test_type)")
end

function hypothesis(io::IO, x::KSOneSampleTest)
  println(io, "Hypothesis to be tested:")
  if x.test_type == :bilateral
    println(io, "  H0:                   F(x) ~ $(x.distribution)")
    println(io, "  H1:                   F(x) ≠ $(x.distribution)")
  elseif x.test_type == :one_sided_right
    println(io, "  H0:                   F(x) ≥ $(x.distribution)")
    println(io, "  H1:                   F(x) < $(x.distribution)")
  elseif x.test_type == :one_sided_left
    println(io, "  H0:                   F(x) ≤ $(x.distribution)")
    println(io, "  H1:                   F(x) > $(x.distribution)")
  end
end

ks_parameters(x::KSOneSampleTest) = (x.α, x.Dn, x.Dcrit)

function ks_critical_one_sample(n::Int64, α::Float64)
  if α == 0.10 && n <= 40
    kolmogorov_crit = [0.95, 0.77639, 0.63604, 0.56522, 0.50945, 0.46799, 0.43607, 0.40962,
      0.38746, 0.36866, 0.35242, 0.33815, 0.32549, 0.31417, 0.30397, 0.29472,
      0.28627, 0.27851, 0.27136, 0.26473, 0.25858, 0.25283, 0.24746, 0.24242,
      0.23768, 0.2332, 0.22898, 0.22497, 0.22117, 0.21756, 0.21412, 0.21085,
      0.20771, 0.20472, 0.20185, 0.1991, 0.19646, 0.19392, 0.19148, 0.18913
    ]
    crit = kolmogorov_crit[n]
  elseif α == 0.05 && n <= 40
    kolmogorov_crit = [0.97500, 0.84189, 0.70760, 0.62394, 0.56328, 0.51926, 0.48342, 0.45427,
      0.43001, 0.40925, 0.39122, 0.37543, 0.36143, 0.34890, 0.33760, 0.32733,
      0.31796, 0.30936, 0.30143, 0.29408, 0.28724, 0.28087, 0.27490, 0.26931,
      0.26404, 0.25907, 0.25438, 0.24993, 0.24571, 0.24170, 0.23788, 0.23424,
      0.23076, 0.22743, 0.22425, 0.22119, 0.21826, 0.21544, 0.21273, 0.21012
    ]
    crit = kolmogorov_crit[n]
  elseif α == 0.01 && n <= 40
    kolmogorov_crit = [0.995, 0.92929, 0.829, 0.73424, 0.66853, 0.61661, 0.57581, 0.54179,
      0.51332, 0.48893, 0.4677, 0.44905, 0.43247, 0.41762, 0.4042, 0.39201,
      0.38086, 0.37062, 0.36117, 0.35241, 0.34427, 0.33666, 0.32954, 0.32286,
      0.31657, 0.31064, 0.30502, 0.29971, 0.29466, 0.28987, 0.2853, 0.28094,
      0.27677, 0.27279, 0.26897, 0.26532, 0.2618, 0.25843, 0.25518, 0.25205
    ]
    crit = kolmogorov_crit[n]
  elseif α == 0.10 && n > 40
    crit = 1.22 / sqrt(n)
  elseif α == 0.05 && n > 40
    crit = 1.3581 / sqrt(n)
  elseif α == 0.01 && n > 40
    crit = 1.6276 / sqrt(n)
  else
    error("α must be one of: 0.10, 0.05, or 0.01, but got $α.")
  end
  return crit
end

function ks_calculated_one_sample(x::AbstractVector{T}, dist::UnivariateDistribution) where {T<:Real}
  n_elements = length(x)
  F_obs = cdf.(dist, sort(x))
  S_x = (1:n_elements) ./ n_elements
  S_x_minus = (0:n_elements-1) ./ n_elements
  Dplus = maximum(S_x .- F_obs)
  Dminus = maximum(F_obs .- S_x_minus)
  Dn = max(Dplus, Dminus)
  return (n_elements, Dplus, Dminus, Dn)
end

function ks_calculated_one_sample(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, dist::UnivariateDistribution)
  n_elements = sum(y)
  F_obs = cumsum(y) ./ n_elements
  S_pdf = cumsum(pdf.(dist, x) * (x[2] - x[1]))
  Dplus = maximum(S_pdf .- F_obs)
  Dminus = maximum(F_obs .- S_pdf)
  Dn = max(Dplus, Dminus)
  return (n_elements, Dplus, Dminus, Dn)
end

function ks_one_sample_acceptance(x::AbstractArray{T}, α::Float64, dist::UnivariateDistribution, test::Symbol=:bilateral) where {T<:Real}
  n_elements, Dplus, Dminus, Dn = ks_calculated_one_sample(x, dist)
  Dcrit = ks_critical_one_sample(n_elements, α)
  return Dn <= Dcrit
end


"""
    ks_one_sample_test(x::AbstractArray{<:Real}, α::Float64, dist::UnivariateDistribution; test::Symbol=:bilateral)

Performs the Kolmogorov–Smirnov one-sample test on the vector `x` with the null hypothesis that the data in `x` 
come from the given distribution `dist`, against the alternative that the sample is not drawn from `dist`, 
for a given significance level α. Valid values for α are 0.10 (90% confidence), 0.05 (95% confidence), or 0.01 (99% confidence).
"""
function ks_one_sample_test(x::AbstractArray{T}, α::Float64, dist::UnivariateDistribution, test::Symbol=:bilateral) where {T<:Real}
  n_elements, Dplus, Dminus, Dn = ks_calculated_one_sample(x, dist)
  Dcrit = ks_critical_one_sample(n_elements, α)
  return KSOneSampleTest(α, n_elements, Dplus, Dminus, Dn, Dcrit, dist, test)
end

"""
    ks_one_sample_test(x::AbstractArray{<:Real}, y::AbstractArray{<:Real}, α::Float64, dist::UnivariateDistribution; test::Symbol=:bilateral)

Performs the Kolmogorov–Smirnov one-sample test on the vector `x` with frequency `y` with the null hypothesis 
that the weighted data in `x` come from the given distribution `dist`, against the alternative that the sample is 
not drawn from `dist`, for a given significance level α. Valid values for α are 0.10 (90% confidence), 0.05 (95% confidence), or 0.01 (99% confidence).
"""
function ks_one_sample_test(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, α::Float64, dist::UnivariateDistribution, test::Symbol=:bilateral)
  n_elements, Dplus, Dminus, Dn = ks_calculated_one_sample(x, y, dist)
  Dcrit = ks_critical_one_sample(n_elements, α)
  return KSOneSampleTest(α, n_elements, Dplus, Dminus, Dn, Dcrit, dist, test)
end

function p_value(x::KSOneSampleTest)
  test_type = x.test_type
  if test_type == :one_sided_left
    return exp(-2 * x.n_elements * x.Dminus^2)
  elseif test_type == :one_sided_right
    return exp(-2 * x.n_elements * x.Dplus^2)
  elseif test_type == :bilateral
    return HypothesisTests.pvalue(Kolmogorov(), sqrt(x.n_elements) * x.Dn; tail=:right)
  else
    throw(ArgumentError("Test=$(test_type) is invalid"))
  end
end

struct KSTwoSampleTest <: KolmogorovSmirnovTest
  α::Float64
  n_elements_x::Int
  n_elements_y::Int
  Dplus::Float64
  Dminus::Float64
  Dn::Float64
  Dcrit::Float64
  test_type::Symbol
end

name_of_test(::KSTwoSampleTest) = "Kolmogorov-Smirnov Two-Sample Test"

function show_params(io::IO, x::KSTwoSampleTest)
  println(io, "  Number of observations:  [$(x.n_elements_x), $(x.n_elements_y)]")
  println(io, "  Test parameters:         $(x.test_type)")
end

function hypothesis(io::IO, x::KSTwoSampleTest)
  println(io, "Hypothesis to be tested:")
  if x.test_type == :bilateral
    println(io, "  H0:                   F(x) = G(x)")
    println(io, "  H1:                   F(x) ≠ G(x)")
  elseif x.test_type == :one_sided_right
    println(io, "  H0:                   F(x) ≥ G(x)")
    println(io, "  H1:                   F(x) < G(x)")
  elseif x.test_type == :one_sided_left
    println(io, "  H0:                   F(x) ≤ G(x)")
    println(io, "  H1:                   F(x) > G(x)")
  end
end

ks_parameters(x::KSTwoSampleTest) = (x.α, x.Dn, x.Dcrit)

function ks_critical_two_sample(n::Int64, α::Float64, test_type::Symbol)
  if α == 0.05 && test_type == :bilateral
    kolmogorov_crit = [NaN, NaN, NaN, 4, 5, 5, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 9, 9, 9,
      10, 10, 10, 10, 10, 11, 11, 11, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13]
    crit = kolmogorov_crit[n]
  elseif α == 0.01 && test_type == :bilateral
    kolmogorov_crit = [NaN, NaN, NaN, NaN, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10, 11,
      11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15]
    crit = kolmogorov_crit[n]
  elseif α == 0.05 && test_type == :one_sided_right
    kolmogorov_crit = [NaN, NaN, 3, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8, 8,
      9, 9, 9, 9, 9, 9, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11]
    crit = kolmogorov_crit[n]
  elseif α == 0.01 && test_type == :one_sided_right
    kolmogorov_crit = [NaN, NaN, NaN, NaN, 5, 6, 6, 6, 7, 7, 8, 8, 8, 8, 9, 9, 9, 10, 10, 10, 10,
      11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14]
    crit = kolmogorov_crit[n]
  elseif test_type == :one_sided_left
    crit = NaN
  else
    error("α must be one of 0.10, 0.05, or 0.01, but got $α.")
  end
  return crit
end

function ks_critical_two_sample(n₁::Int64, n₂::Int64, α::Float64, test_type::Symbol)
  if α == 0.05 && test_type == :bilateral
    return 1.36 * sqrt((n₁ + n₂) / (n₁ * n₂))
  elseif α == 0.01 && test_type == :bilateral
    return 1.63 * sqrt((n₁ + n₂) / (n₁ * n₂))
  elseif α == 0.05 && test_type == :one_sided_right
    return 5.991
  elseif α == 0.01 && test_type == :one_sided_right
    return 9.210
  elseif α == 0.05 && test_type == :one_sided_left
    return 5.991
  elseif α == 0.01 && test_type == :one_sided_left
    return 9.210
  else
    error("α must be one of 0.10 or 0.05, but got $α.")
  end
end

function ks_calculated_two_sample(x::AbstractVector{T}, y::AbstractVector{T}, test_type::Symbol) where {T<:Real}
  n_elements_x, n_elements_y = length(x), length(y)
  permutation = sortperm([x; y])
  weights = [ones(n_elements_x) / n_elements_x; -ones(n_elements_y) / n_elements_y][permutation]
  cumulative = cumsum(weights)
  Dplus = maximum(cumulative)
  Dminus = -minimum(cumulative)
  Dn = max(Dplus, Dminus)
  if test_type == :one_sided_right && max(n_elements_x, n_elements_y) <= 40
    return Dplus, Dminus, Dplus
  elseif test_type == :one_sided_left && max(n_elements_x, n_elements_y) <= 40
    return Dplus, Dminus, Dminus
  elseif test_type == :one_sided_right && max(n_elements_x, n_elements_y) > 40
    Dn = 4 * (Dplus)^2 * ((n_elements_x * n_elements_y) / (n_elements_x + n_elements_y))
    return Dplus, Dminus, Dn
  elseif test_type == :one_sided_left && max(n_elements_x, n_elements_y) > 40
    Dn = 4 * (Dminus)^2 * ((n_elements_x * n_elements_y) / (n_elements_x + n_elements_y))
    return Dplus, Dminus, Dn
  else
    return Dplus, Dminus, Dn
  end
end

function ks_two_sample_test(x::AbstractVector{T}, y::AbstractVector{T}, α::Float64=0.05, test_type::Symbol=:bilateral) where {T<:Real}
  n_elements_x, n_elements_y = length(x), length(y)
  Dplus, Dminus, Dn = ks_calculated_two_sample(x, y, test_type)
  if max(n_elements_x, n_elements_y) <= 40
    Dcrit = ks_critical_two_sample(n_elements_x, α, test_type) / ((n_elements_x + n_elements_y) / 2)
  else
    Dcrit = ks_critical_two_sample(n_elements_x, n_elements_y, α, test_type)
  end
  return KSTwoSampleTest(α, n_elements_x, n_elements_y, Dplus, Dminus, Dn, Dcrit, test_type)
end

function p_value(x::KSTwoSampleTest)
  n = (x.n_elements_x * x.n_elements_y) / (x.n_elements_x + x.n_elements_y)
  test_type = x.test_type
  if test_type == :one_sided_left
    return exp(-2 * n * x.Dminus^2)
  elseif test_type == :one_sided_right
    return exp(-2 * n * x.Dplus^2)
  elseif test_type == :bilateral
    return HypothesisTests.pvalue(Kolmogorov(), sqrt(n) * x.Dn; tail=:right)
  else
    throw(ArgumentError("Test=$(test_type) is invalid"))
  end
end


