function _nulldeviance(x::Vector{<:Real})
  mean_x = mean(x)
  out = Vector{Float64}(undef, length(x))
  @simd for i in eachindex(x)
    out[i] = abs2(x[i] - mean_x)
  end
  return sum(out)
end

_p_result(test::HypothesisTests.HypothesisTest) = pvalue(test) > 0.05 ? true : false

"""
Function to calculate various statistical criteria for evaluating regression models.

This function computes several metrics that are useful for comparing the performance of 
different regression models, including goodness-of-fit, residual error metrics, and model quality measures.

Arguments:
- `model::FittedLinearModel`: A fitted linear regression model.

Returns:
- `parameters::Matrix{Float64}`: A matrix containing the calculated metrics.
"""
function _criteria_parameters(model::FittedLinearModel)::Matrix{Float64}
  # Number of observations (n) in the dataset
  n = model.data[1] |> length
  # Model coefficients
  coefficients = model.β
  ncoef = length(coefficients)
  # Degrees of freedom for residuals (dof_resid)
  dof_resid = n - ncoef
  # Actual observed values (response variable)
  y = model.data[1]
  # Predicted values from the model
  ŷ = predict(model)
  # Residuals: the difference between observed and predicted values
  residual = y - ŷ
  # Deviance: sum of squared residuals (SSR)
  SSR = residual ⋅ residual
  # Mean Squared Error (MSE)
  MSE = SSR / n
  # Root Mean Squared Error (RMSE)
  RMSE = √MSE
  # Mean Absolute Error (MAE): average absolute residual
  MAE = mean(abs.(residual))
  # Standard error of the estimate (Syx) as a percentage of the mean response
  syx_in_percentage = √(SSR / dof_resid) / mean(y) * 100
  # Coefficient of determination (R²): proportion of variance explained
  r_2 = 1 - SSR / _nulldeviance(y)
  # Adjusted R²: adjusted for the number of predictors
  adj_r² = 1 - (1 - r_2) * (n - 1) / dof_resid
  # Log-likelihood of the model
  loglike = -n / 2 * (log(2π * MSE) + 1)
  # Akaike Information Criterion (AIC): model quality measure
  AIC = -2 * loglike + 2 * ncoef
  # Bayesian Information Criterion (BIC): penalizes model complexity more than AIC
  BIC = -2 * loglike + log(n) * ncoef

  # Normality test for residuals using Kolmogorov-Smirnov test
  normality_p = try
    if n < 100
      ExactOneSampleKSTest(HypothesisTests.ksstats(residual, fit_mle(Normal, residual))...) |> _p_result
    else
      ApproximateOneSampleKSTest(HypothesisTests.ksstats(residual, fit_mle(Normal, residual))...) |> _p_result
    end
  catch
    false
  end

  # Test for coefficient significance using p-values
  dispersion = rmul!(inv(model.chol), model.σ²)
  standard_errors = sqrt.(diag(dispersion))
  t_values = coefficients ./ standard_errors
  p_values = ccdf.(Ref(FDist(1, dof_resid)), abs2.(t_values))

  # Check if all coefficients are significant at the 0.05 level
  coefs_significant = all(p_values .< 0.05) ? true : false

  # Return a matrix containing the calculated parameters
  parameters = [r_2 adj_r² syx_in_percentage RMSE MAE AIC BIC normality_p coefs_significant]
  return parameters
end

function _calculate_ranks(ct::DataFrame, selected_criteria::Vector{Symbol})

  n = size(ct, 1)
  ranks = Dict()

  # Calculate ranks for each criterion
  if :r2 in selected_criteria
    ranks[:r2] = competerank(ct[!, "r2"], rev=true)
  end
  if :adjr2 in selected_criteria
    ranks[:adjr2] = competerank(ct[!, "adjr2"], rev=true)
  end
  if :pvalue in selected_criteria
    ranks[:pvalue] = competerank(ct[!, "pvalue"], rev=true)
  end
  if :syx in selected_criteria
    ranks[:syx] = competerank(ct[!, "syx"])
  end
  if :rmse in selected_criteria
    ranks[:rmse] = competerank(ct[!, "rmse"])
  end
  if :mae in selected_criteria
    ranks[:mae] = competerank(ct[!, "mae"])
  end
  if :aic in selected_criteria
    ranks[:aic] = competerank(ct[!, "aic"])
  end
  if :bic in selected_criteria
    ranks[:bic] = competerank(ct[!, "bic"])
  end
  if :normality in selected_criteria
    # Penalize non-normal models with a higher rank
    normality_ranks = competerank(ct[!, "normality"], rev=true)
    penalized_non_normal_ranks = [normality_ranks[i] == 1 ? 1 : normality_ranks[i] * n for i in 1:n]
    ranks[:normality] = penalized_non_normal_ranks
  end
  if :significance in selected_criteria
    # Penalize non-significance models with a higher rank
    significance_ranks = competerank(ct[!, "significance"], rev=true)
    penalized_non_significance_ranks = [significance_ranks[i] == 1 ? 1 : significance_ranks[i] * n for i in 1:n]
    ranks[:significance] = penalized_non_significance_ranks
  end

  # Combine ranks into a single score
  combined_rank = sum([ranks[crit] for crit in selected_criteria])

  return combined_rank
end

"""
    criteria_table(model::Vector{<:FittedLinearModel}, criteria::Symbol...; best::Union{Bool,Int}=10)

The `criteria_table` function evaluates and ranks multiple regression models based on specified criteria. 
  It generates a comprehensive table of performance metrics for each model, calculates ranks for these 
    metrics, and combines them into a final score. The function allows for flexible selection of 
    evaluation criteria and can return either all models or only the top models based on the combined ranking.

# Parameters:
- `model`: 
  The regression model(s) to be evaluated and compared. This parameter can accept:
  - **Single Linear Regression Model (`FittedLinearModel`)**:
    Evaluates a single linear regression model.
  - **Vector of Linear Regression Models (`Vector{<:FittedLinearModel}`)**:
    Evaluates and compares multiple linear regression models.
    
- `criteria::Symbol...`: 
  A variable number of symbols representing the evaluation criteria to include. Possible values include:
  - `:r2`: R², a measure of the model's explanatory power, representing the proportion of variance in the dependent variable explained by the predictors.
  - `:adjr2`: Adjusted R², a measure of the model's explanatory power, adjusted for the number of predictors.
  - `:syx`: Standard error of the estimate (Syx%) expressed as a percentage of the mean of the dependent 
  variable (y), indicating the precision of the model's predictions.
  - `:rmse`: Root Mean Squared Error, indicating the average magnitude of residuals.
  - `:mae`: Mean Absolute Error, another accuracy measure based on average absolute residuals.
  - `:aic`: Akaike Information Criterion, balancing goodness of fit with model complexity.
  - `:bic`: Bayesian Information Criterion, penalizes model complexity more than AIC;
  - `:normality`: Assesses the normality of residuals using the Kolmogorov-Smirnov test, ensuring that residuals follow a normal distribution.
  - `:significance`: Evaluates whether model coefficients are statistically significant.
  - `:all`: Selects all available criteria for evaluation.

  If `:all` is specified, the function will use every available criterion: [:r2, :adjr2, :syx, :rmse, :mae, :aic, :bic, :normality, :significance]. 
  If no criteria are specified, the function defaults to the recommended criteria: [:adjr2, :syx, :aic, :bic, :normality, :significance].

- `best::Union{Bool, Int}=10`: 
  Specifies the number of top models to return based on the combined ranking.
  - `false`: Returns the full table with all models ranked.
  - Integer value: If less than the total number of models, returns only the top `best` models.

# Returns:
- `DataFrame`: 
    A sorted table with the evaluated models and their respective metrics. Includes a combined rank 
    based on the selected criteria.

# Examples:
- **Single Model**: 
  `criteria_table(model, :adjr2, :bic)`
  
- **Vector of Models**:
  `criteria_table([model1, model2], :aic, :mae)`
"""
function criteria_table(
  model::Vector{<:FittedLinearModel},
  criteria::Symbol...;
  best::Union{Bool,Int}=10
)

  # Define allowed fields for criteria
  allowed_fields = [:r2, :adjr2, :syx, :rmse, :mae, :aic, :bic, :normality, :significance]

  # Determine selected criteria
  if isempty(criteria)
    # Use default criteria if none are specified
    selected_criteria = [:adjr2, :syx, :aic, :bic, :normality, :significance]
  elseif :all in criteria
    # If :all is included, use all fields except :all itself
    selected_criteria = setdiff(allowed_fields, [:all])
  else
    # Use explicitly specified criteria
    selected_criteria = Symbol[s for s in criteria]
  end

  # Validate selected criteria
  if !issubset(Set(selected_criteria), Set(setdiff(allowed_fields, [:all])))
    not_allowed = join(setdiff(selected_criteria, allowed_fields), ", :")
    allowed_msg = "\nAllowed fields are: :" * join(allowed_fields, ", :")
    throw(ArgumentError(":$not_allowed not allowed." * allowed_msg))
  end

  # Generate criteria parameters for each model
  criteria_params = vcat(_criteria_parameters.(model)...)

  # Create a DataFrame from the full criteria parameters
  ct = DataFrame(criteria_params, allowed_fields)

  # Filter the DataFrame columns based on the selected criteria
  ct = select(ct, selected_criteria)

  # Insert the model objects into the DataFrame
  insertcols!(ct, 1, "model" => model)

  # Combine ranks into a single score
  combined_rank = _calculate_ranks(ct, selected_criteria)

  # Insert the combined rank into the DataFrame
  insertcols!(ct, 2, :rank => combined_rank)

  # Sort the DataFrame by the combined rank
  sort!(ct, :rank)

  # If 'best' is false, return the full DataFrame
  if best === false
    return ct
  elseif best < length(model)
    # If 'best' is less than the number of models, return the top 'best' models
    top_models = ct[1:best, 1]
    return criteria_table(top_models, criteria...; best=false) # Re-run with selected models
  else
    # Otherwise, return the full DataFrame
    return ct
  end
end

criteria_table(model::FittedLinearModel, criteria::Symbol...) = criteria_table([model], criteria...)

"""
    criteria_selection(model::Vector{<:FittedLinearModel}, criteria::Symbol...)

The `criteria_selection` function evaluates and ranks a vector of regression models based on specified 
  criteria, returning the best model according to the combined ranking.

# Parameters:
- `model::Vector{<:FittedLinearModel}`: 
  A vector of linear regression models to be evaluated and compared.

- `criteria::Symbol...`: 
  A variable number of symbols representing the evaluation criteria to include. Possible values include:
  - `:adjr2`: Adjusted R², a measure of the model's explanatory power, adjusted for the number of predictors.
  - `:syx`: Standard error of the estimate as a percentage of the mean of `y`.
  - `:rmse`: Root Mean Squared Error, indicating the average magnitude of residuals.
  - `:mae`: Mean Absolute Error, another accuracy measure based on average absolute residuals.
  - `:aic`: Akaike Information Criterion, balancing goodness of fit with model complexity.
  - `:significance`: Evaluates whether model coefficients are statistically significant.
  - `:normality`: Assesses the normality of residuals, an assumption in linear regression.
  - `:homoscedasticity`: Checks for constant variance in residuals, another key regression assumption.
  
  If no criteria are specified, the function will use all available criteria by default.

# Returns:
- `FittedLinearModel`: 
  The best model based on the combined ranking of the specified criteria.
"""
criteria_selection(model::Vector{<:FittedLinearModel}, criteria::Symbol...) = criteria_table(model, criteria..., best=5)[1, 1]