"""
Predict values using a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.
- `data::AbstractDataFrame`: The data frame containing the new data for prediction.

# Returns
- `Vector{Float64}`: The predicted values.
"""
function predict(fitted_model::FittedLinearModel, data::AbstractDataFrame)
  x, nonmissings = StatsModels.missing_omit(columntable(data), fitted_model.formula.rhs)
  X = modelmatrix(fitted_model.formula.rhs, x)
  ŷ = X * fitted_model.β
  if isa(fitted_model.formula.lhs, FunctionTerm)
    ŷ = _predict(ŷ, x[1], fitted_model.σ², nameof(fitted_model.formula.lhs.f))
  end
  length(unique(nonmissings)) == 1 ? ŷ : StatsModels._return_predictions(Tables.materializer(data), ŷ, nonmissings, length(nonmissings))
end

"""
Predict a single value using a fitted linear model for a single row of data.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.
- `data::DataFrameRow`: A single row of data for prediction.

# Returns
- `Float64`: The predicted value.
"""
function predict(fitted_model::FittedLinearModel, data::DataFrameRow)
  x, nonmissings = StatsModels.missing_omit(columntable(DataFrame(data)), fitted_model.formula.rhs)
  if isempty(x[1])
    return missing
  else
    X = modelmatrix(fitted_model.formula.rhs, x)
    ŷ = X * fitted_model.β
    if isa(fitted_model.formula.lhs, FunctionTerm)
      ŷ = _predict(ŷ, x[1], fitted_model.σ², nameof(fitted_model.formula.lhs.f))
    end
    return ŷ[1]
  end
end

"""
Predict values using the data from the fitted model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Vector{Float64}`: The predicted values.
"""
@inline predict(fitted_model::FittedLinearModel) = predict(fitted_model, fitted_model.data)

"""
Calculate residuals of a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Vector{Float64}`: The residuals.
"""
residuals(fitted_model::FittedLinearModel) = fitted_model.data[!, 1] - predict(fitted_model)

"""
Get coefficients of a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Vector{Float64}`: The coefficients.
"""
coef(fitted_model::FittedLinearModel) = fitted_model.β

"""
Get the number of coefficients of a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Int`: The number of coefficients.
"""
n_coef(fitted_model::FittedLinearModel) = length(coef(fitted_model))

"""
Get names of the coefficients of a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Vector{String}`: The names of the coefficients.
"""
@inline coefnames(fitted_model::FittedLinearModel) = vcat("β0", string.(StatsModels.coefnames(fitted_model.formula.rhs.terms[2:end])))

"""
Get degrees of freedom of a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Int`: The degrees of freedom.
"""
dof(fitted_model::FittedLinearModel) = n_coef(fitted_model) + 1

"""
Get the number of observations used in fitting the model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Int`: The number of observations.
"""
nobs(fitted_model::FittedLinearModel) = size(fitted_model.data, 1)

"""
Get degrees of freedom for residuals of a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Int`: The degrees of freedom for residuals.
"""
dof_residual(fitted_model::FittedLinearModel) = nobs(fitted_model) - dof(fitted_model) + 1

"""
Get the model matrix of a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Matrix{Float64}`: The model matrix.
"""
StatsModels.modelmatrix(fitted_model::FittedLinearModel) = modelmatrix(fitted_model.formula, fitted_model.data)

"""
Calculate the dispersion of a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Matrix{Float64}`: The dispersion matrix.
"""
dispersion(fitted_model::FittedLinearModel) = rmul!(inv(fitted_model.chol), fitted_model.σ²)

"""
Calculate the standard error of the coefficients of a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Vector{Float64}`: The standard errors.
"""
stderror(fitted_model::FittedLinearModel) = sqrt.(diag(dispersion(fitted_model)))

"""
Create a table of coefficients for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `DataFrame`: The coefficient table.
"""
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

"""
Calculate the null deviance for a vector of values.

# Arguments
- `x::Vector{<:Real}`: The vector of values.

# Returns
- `Float64`: The null deviance.
"""
# function StatsBase.nulldeviance(x::Vector{<:Real})
#   out = similar(x)
#   @inbounds for i in eachindex(x)
#     out[i] = abs2(x[i] - mean(x))
#   end
#   return sum(out)
# end
function StatsBase.nulldeviance(x::Vector{<:Real})
  mean_x = mean(x)
  out = Vector{Float64}(undef, length(x))
  @inbounds for i in eachindex(x)
    out[i] = abs2(x[i] - mean_x)
  end
  return sum(out)
end

"""
Calculate the deviance for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The deviance.
"""
function StatsBase.deviance(fitted_model::FittedLinearModel)
  resid = residuals(fitted_model)
  return resid ⋅ resid
end

"""
Calculate the null deviance for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The null deviance.
"""
StatsBase.nulldeviance(fitted_model::FittedLinearModel) = nulldeviance(fitted_model.data[!, 1])

"""
Calculate the R-squared value for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The R-squared value.
"""
r2(fitted_model::FittedLinearModel) = 1 - deviance(fitted_model) / nulldeviance(fitted_model)

"""
Calculate the adjusted R-squared value for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The adjusted R-squared value.
"""
function adjr2(fitted_model::FittedLinearModel)
  n = nobs(fitted_model)
  p = n_coef(fitted_model)
  1 - (1 - r2(fitted_model)) * (n - 1) / (n - p)
end

"""
Calculate the log-likelihood for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The log-likelihood.
"""
loglikelihood(fitted_model::FittedLinearModel) = -nobs(fitted_model) / 2 * (log(2π * deviance(fitted_model) / nobs(fitted_model)) + 1)

"""
Calculate the Akaike Information Criterion (AIC) for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The AIC value.
"""
aic(fitted_model::FittedLinearModel) = -2 * loglikelihood(fitted_model) + 2 * dof(fitted_model)

"""
Calculate the corrected Akaike Information Criterion (AICc) for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The AICc value.
"""
function aicc(fitted_model::FittedLinearModel)
  k = dof(fitted_model)
  n = nobs(fitted_model)
  -2 * loglikelihood(fitted_model) + 2 * k + 2 * k * (k + 1) / (n - k - 1)
end

"""
Calculate the Bayesian Information Criterion (BIC) for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The BIC value.
"""
bic(fitted_model::FittedLinearModel) = -2 * loglikelihood(fitted_model) + dof(fitted_model) * log(nobs(fitted_model))

"""
Calculate the standard error of the estimate for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The standard error of the estimate.
"""
syx(fitted_model::FittedLinearModel) =  √(deviance(fitted_model) / dof_residual(fitted_model))

"""
Calculate the standard error of the estimate as a percentage for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Float64`: The standard error of the estimate as a percentage.
"""
syx_in_percentage(fitted_model::FittedLinearModel) = syx(fitted_model) / mean(fitted_model.data[!, 1]) * 100

"""
Test for normality of residuals in a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `ExactOneSampleKSTest`: The result of the Kolmogorov-Smirnov test.
"""
function normality(fitted_model::FittedLinearModel)
  residual = residuals(fitted_model)
  ExactOneSampleKSTest(ksstats(residual, fit_mle(Normal, residual))...)
end

"""
Test for homoscedasticity of residuals in a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `WhiteTest`: The result of the White test.
"""
@inline homoscedasticity(fitted_model::FittedLinearModel) = WhiteTest([ones(nobs(fitted_model)) predict(fitted_model)], residuals(fitted_model))

"""
Check the result of a hypothesis test.

# Arguments
- `test::HypothesisTests.HypothesisTest`: The hypothesis test result.

# Returns
- `Bool`: True if the p-value is greater than 0.05, otherwise false.
"""
@inline p_result(test::HypothesisTests.HypothesisTest) = pvalue(test) > 0.05 ? true : false

"""
Calculate various criteria parameters for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Matrix{Float64}`: The matrix of criteria parameters.
"""
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
  normal = try
      ExactOneSampleKSTest(HypothesisTests.ksstats(residual, fit_mle(Normal, residual))...) |> p_result
    catch 
      false
  end
  homosced = try
    WhiteTest([ones(nobs) ŷ], residual) |> p_result
  catch 
    false
  end
  parameters = [round(fitted_model.RMSE, digits = 4) syx_in_percentage adjr2 normal homosced]
  return parameters
end

"""
Create a criteria table for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `DataFrame`: The criteria table.
"""
function criteria_table(fitted_model::FittedLinearModel) :: DataFrame
  ct = DataFrame(_criteria_parameters(fitted_model), ["RMSE", "Syx%", "Adj. R²", "Normality", "Homoscedasticity"])
  insertcols!(ct, 1, "Adj. Model" => fitted_model)
  return ct
end

"""
Create a criteria table for multiple fitted linear models.

# Arguments
- `fitted_model::Vector{<:FittedLinearModel}`: The vector of fitted linear models.

# Returns
- `DataFrame`: The criteria table.
"""
# function criteria_table(fitted_model::Vector{<:FittedLinearModel}) :: DataFrame
#   ct = DataFrame(vcat(_criteria_parameters.(fitted_model)...), ["RMSE", "Syx%", "Adj. R²", "Normality", "Homoscedasticity"])
#   insertcols!(ct, 1, "Adj. Model" => fitted_model)
#   insertcols!(ct, 2, :Rank => .+(competerank(ct[!, 2]), map(i -> competerank(ct[!, i], rev = true), 4:6)...))
#   sort!(ct, :Rank)
#   return ct
# end
function criteria_table(fitted_model::Vector{<:FittedLinearModel}) :: DataFrame
  # Generate the criteria parameters for each model
  criteria_params = _criteria_parameters.(fitted_model)
  
  # Create a DataFrame from the criteria parameters
  ct = DataFrame(vcat(criteria_params...), ["RMSE", "Syx%", "Adj. R²", "Normality", "Homoscedasticity"])
  
  # Insert the models and calculate the ranks
  insertcols!(ct, 1, "Adj. Model" => fitted_model)
  
  # Calculate the ranks for each criterion
  rmse_rank = competerank(ct[!, "RMSE"])
  syx_rank = competerank(ct[!, "Syx%"])
  adj_r2_rank = competerank(ct[!, "Adj. R²"], rev = true)
  normality_rank = competerank(ct[!, "Normality"], rev = true)
  homoscedasticity_rank = competerank(ct[!, "Homoscedasticity"], rev = true)
  
  # Combine ranks with weighted priorities (prioritize normality and homoscedasticity)
  combined_rank = rmse_rank .+ syx_rank .+ adj_r2_rank .+ 10 * normality_rank .+ 10 * homoscedasticity_rank
  
  # Insert the combined rank into the DataFrame
  insertcols!(ct, 2, :Rank => combined_rank)
  
  # Sort the DataFrame by rank
  sort!(ct, :Rank)
  
  return ct
end






"""
Calculate the confidence intervals for a fitted linear model.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted linear model.

# Returns
- `Vector{Float64}`: The confidence intervals.
"""
function confidence_interval(fitted_model::FittedLinearModel)
  model_matrix = modelmatrix(fitted_model)
  R = fitted_model.chol.U
  variance = deviance(fitted_model) / dof_residual(fitted_model)
  residvar = ones(size(model_matrix, 2)) * variance
  tvalue = -quantile(TDist(dof_residual(fitted_model)), 0.025)
  confidence = tvalue * .√((model_matrix / R).^2 * residvar)
  return confidence
end