function predict(fitted_model::FittedLinearModel, data::AbstractDataFrame)
  x, nonmissings = StatsModels.missing_omit(columntable(data), fitted_model.formula.rhs)
  X = modelmatrix(fitted_model.formula.rhs, x)
  ŷ = X * fitted_model.β
  if isa(fitted_model.formula.lhs, FunctionTerm)
    ŷ = _predict(ŷ, x[1], fitted_model.σ², nameof(fitted_model.formula.lhs.f))
  end
  length(unique(nonmissings)) == 1 ? ŷ : _return_predictions(Tables.materializer(data), ŷ, nonmissings, length(nonmissings))
end

function predict(fitted_model::FittedLinearModel, data::DataFrameRow)
  x, nonmissings = StatsModels.missing_omit(columntable(DataFrame(data)), fitted_model.formula.rhs)
  X = modelmatrix(fitted_model.formula.rhs, x)
  ŷ = X * fitted_model.β
  if isa(fitted_model.formula.lhs, FunctionTerm)
    ŷ = _predict(ŷ, x[1], fitted_model.σ², nameof(fitted_model.formula.lhs.f))
  end
  return ŷ[1]
end

@inline predict(fitted_model::FittedLinearModel) = predict(fitted_model, fitted_model.data)

residuals(fitted_model::FittedLinearModel) = fitted_model.data[!, 1] - predict(fitted_model)

coef(fitted_model::FittedLinearModel) = fitted_model.β

n_coef(fitted_model::FittedLinearModel) = length(coef(fitted_model))

@inline coefnames(fitted_model::FittedLinearModel) = vcat("β0", string.(StatsModels.coefnames(fitted_model.formula.rhs.terms[2:end])))

dof(fitted_model::FittedLinearModel) = n_coef(fitted_model) + 1

nobs(fitted_model::FittedLinearModel) = size(fitted_model.data, 1)

dof_residual(fitted_model::FittedLinearModel) = nobs(fitted_model) - dof(fitted_model) + 1

StatsModels.modelmatrix(fitted_model::FittedLinearModel) = modelmatrix(fitted_model.formula, fitted_model.data)

dispersion(fitted_model::FittedLinearModel) = rmul!(inv(fitted_model.chol), fitted_model.σ²)

stderror(fitted_model::FittedLinearModel) = sqrt.(diag(dispersion(fitted_model)))

function coef_table(fitted_model::FittedLinearModel) :: DataFrame
  β = coef(fitted_model)
  std_error = stderror(fitted_model)
  t_value = β ./ std_error
  p_value = ccdf.(Ref(FDist(1, dof_residual(fitted_model))), abs2.(t_value))
  limit_interval = std_error * quantile(TDist(dof_residual(fitted_model)), 0.05 / 2)
  lower_limit = β + limit_interval
  upper_limit = β - limit_interval
  DataFrame(
    "Coef. Names" => coefnames(fitted_model), "Coef. Values" => β,
    "Std. Error" =>  std_error, "T Value" =>  t_value, "Pr(>|t|)" =>  p_value,
    "Lower 95%" =>  lower_limit,"Upper 95%" =>  upper_limit
  )
end

function StatsBase.nulldeviance(x::Vector{<:Real})
  out = similar(x)
  @inbounds for i in eachindex(x)
    out[i] = abs2(x[i] - mean(x))
  end
  return sum(out)
end

function StatsBase.deviance(fitted_model::FittedLinearModel)
  resid = residuals(fitted_model)
  return resid ⋅ resid
end

StatsBase.nulldeviance(fitted_model::FittedLinearModel) = nulldeviance(fitted_model.data[!, 1])

r2(fitted_model::FittedLinearModel) = 1 - deviance(fitted_model) / nulldeviance(fitted_model)

function adjr2(fitted_model::FittedLinearModel)
  n = nobs(fitted_model)
  p = n_coef(fitted_model)
  1 - (1 - r2(fitted_model)) * (n - 1) / (n - p)
end

loglikelihood(fitted_model::FittedLinearModel) = -nobs(fitted_model) / 2 * (log(2π * deviance(fitted_model) / nobs(fitted_model)) + 1)

aic(fitted_model::FittedLinearModel) = -2 * loglikelihood(fitted_model) + 2 * dof(fitted_model)
# Akaike information criterion not biased for small samples, when (n/p) < 40
function aicc(fitted_model::FittedLinearModel)
  k = dof(fitted_model)
  n = nobs(fitted_model)
  -2 * loglikelihood(fitted_model) + 2 * k + 2 * k * (k + 1) / (n - k - 1)
end

bic(fitted_model::FittedLinearModel) = -2 * loglikelihood(fitted_model) + dof(fitted_model) * log(nobs(fitted_model))

"The standard error of the estimate is a measure of the accuracy of predictions."
syx(fitted_model::FittedLinearModel) =  √(deviance(fitted_model) / dof_residual(fitted_model))

syx_in_percentage(fitted_model::FittedLinearModel) = syx(fitted_model) / mean(fitted_model.data[!, 1]) * 100

function normality(fitted_model::FittedLinearModel)
  residual = residuals(fitted_model)
  ExactOneSampleKSTest(ksstats(residual, fit_mle(Normal, residual))...)
end

@inline homoscedasticity(fitted_model::FittedLinearModel) = WhiteTest([ones(nobs(fitted_model)) predict(fitted_model)], residuals(fitted_model))

@inline p_result(test::HypothesisTests.HypothesisTest) = pvalue(test) > 0.05 ? true : false

function _criteria_parameters(fitted_model::FittedLinearModel) :: Matrix{Float64}
  nobs = size(fitted_model.data, 1)
  n_coef = length(fitted_model.β)
  dof_residual = nobs - (n_coef + 1) + 1
  ŷ = predict(fitted_model)
  residual = fitted_model.data[!, 1] - ŷ
  deviance = residual ⋅ residual
  syx_in_percentage = round(√(deviance / dof_residual) / mean(fitted_model.data[!, 1]) * 100, digits = 2)
  r2 = 1 - deviance / nulldeviance(fitted_model.data[!, 1])
  adjr2 = round(1 - (1 - r2) * (nobs - 1) / (nobs - n_coef), digits = 4)
  normal = ExactOneSampleKSTest(HypothesisTests.ksstats(residual, fit_mle(Normal, residual))...) |> p_result
  homosced = WhiteTest([ones(nobs) ŷ], residual) |> p_result
  parameters = [round(fitted_model.RMSE, digits = 4) syx_in_percentage adjr2 normal homosced]
  return parameters
end

function criteria_table(fitted_model::FittedLinearModel) :: DataFrame
  ct = DataFrame(_criteria_parameters(fitted_model), ["RMSE", "Syx%", "Adj. R²", "Normality", "Homoscedasticity"])
  insertcols!(ct, 1, "Adj. Model" => fitted_model)
  return ct
end

function criteria_table(fitted_model::Vector{<:FittedLinearModel}) :: DataFrame
  ct = DataFrame(vcat(_criteria_parameters.(fitted_model)...), ["RMSE", "Syx%", "Adj. R²", "Normality", "Homoscedasticity"])
  insertcols!(ct, 1, "Adj. Model" => fitted_model)
  insertcols!(ct, 2, :Rank => .+(competerank(ct[!, 2]), map(i -> competerank(ct[!, i], rev = true), 4:6)...))
  sort!(ct, :Rank)
  return ct
end

function confidence_interval(fitted_model::FittedLinearModel)
  model_matrix = modelmatrix(fitted_model)
  R = fitted_model.chol.U
  variance = deviance(fitted_model) / dof_residual(fitted_model)
  residvar = ones(size(model_matrix, 2)) * variance
  tvalue = -quantile(TDist(dof_residual(fitted_model)), 0.025)
  confidence = tvalue * .√((model_matrix / R).^2 * residvar)
  return confidence
end