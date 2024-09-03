function StatsBase.nulldeviance(x::Vector{<:Real})
  mean_x = mean(x)
  out = Vector{Float64}(undef, length(x))
  @inbounds @simd for i in eachindex(x)
    out[i] = abs2(x[i] - mean_x)
  end
  return sum(out)
end

@inline p_result(test::HypothesisTests.HypothesisTest) = pvalue(test) > 0.05 ? true : false

# Function to calculate various statistical criteria for evaluating regression models
function _criteria_parameters(model::TableRegressionModel{<:LinearModel}) :: Matrix{Float64}
  # Number of observations in the model
  n = nobs(model)
  # Degrees of freedom for residuals
  dof_resid = dof_residual(model)
  # The actual observed values (response variable)
  y = model.mf.data[1]
  # Predicted values from the model
  ŷ = prediction(model)
  # Residuals: the difference between observed and predicted values
  residual = y - ŷ
  # Deviance: sum of squared residuals
  devi = residual ⋅ residual
  # Calculate the Root Mean Squared Error (RMSE)
  RMSE = √(devi / n)
  # Calculate the Mean Absolute Error (MAE)
  MAE = mean(abs.(residual))
  # Standard error of the estimate (Syx) expressed as a percentage of the mean of y
  syx_in_percentage = √(devi / dof_resid) / mean(y) * 100
  # Coefficient of determination (R²)
  r_2 = 1 - devi / nulldeviance(y)
  # Adjusted R²: adjusted for the number of predictors in the model
  adj_r² = 1 - (1 - r_2) * (n - 1) / dof_resid
  # Log-likelihood of the model
  loglike = -n / 2 * (log(2π * devi / n) + 1)
  # Akaike Information Criterion (AIC): a measure of model quality
  AIC = -2 * loglike + 2 * dof_resid
  # Test for coefficient significance
  cc = coef(model)
  se = stderror(model)
  tt = cc ./ se
  p_values = ccdf.(Ref(FDist(1, dof_resid)), abs2.(tt))
  # Check if all coefficients are significant at the 0.05 level
  coefs_significant = all(p_values .< 0.05) ? true : false
  # Normality test for residuals using the Kolmogorov-Smirnov test
  normal = try
    if n < 100
      ExactOneSampleKSTest(HypothesisTests.ksstats(residual, fit_mle(Normal, residual))...) |> p_result
    else
      ApproximateOneSampleKSTest(HypothesisTests.ksstats(residual, fit_mle(Normal, residual))...) |> p_result
    end
    catch
      false
  end
  
  # Homoscedasticity test using the White test
  homoscedasticity = try
    WhiteTest(modelmatrix(model), residual) |> p_result
  catch 
    false
  end
  # Return a matrix containing the calculated parameters
  parameters = [adj_r² syx_in_percentage RMSE MAE AIC coefs_significant normal homoscedasticity]

  return parameters
end

# Function to calculate various statistical criteria for evaluating regression models
function _criteria_parameters(models::GroupedLinearModel)
  # The actual observed values (response variable)
  y = models.group_data[!, 1]
  # grouped_models values
  grouped_models = models.grouped_models |> values |> collect
  # Number of observations in the models
  n = length(y)
  # coef
  cc = coef.(grouped_models)
  n_coefs = sum(length.(cc))
  # Degrees of freedom for residuals
  dof_resid = n - n_coefs
  # Predicted values from the models
  ŷ = prediction(models) |> skipmissing |> collect
  # Residuals: the difference between observed and predicted values
  residual = y - ŷ
  # # Deviance: sum of squared residuals
  devi = residual ⋅ residual
  # # Calculate the Root Mean Squared Error (RMSE)
  RMSE = √(devi / n)
  # # Calculate the Mean Absolute Error (MAE)
  MAE = mean(abs.(residual))
  # # Standard error of the estimate (Syx) expressed as a percentage of the mean of y
  syx_in_percentage = √(devi / dof_resid) / mean(y) * 100
  # # Coefficient of determination (R²)
  r_2 = 1 - devi / nulldeviance(y)
  # # Adjusted R²: adjusted for the number of predictors in the models
  adj_r² = 1 - (1 - r_2) * (n - 1) / dof_resid
  # # Log-likelihood of the models
  loglike = -n / 2 * (log(2π * devi / n) + 1)
  # Akaike Information Criterion (AIC): a measure of models quality
  AIC = -2 * loglike + 2 * dof_resid
  # Test for coefficient significance
  se = stderror.(grouped_models)
  p_values = [ccdf.(Ref(FDist(1, dof_residual(m))), abs2.(cc[i] ./ se[i])) for (i, m) in enumerate(grouped_models)]
  # # Check if all coefficients are significant at the 0.05 level
  coefs_significant = all(vcat(p_values...) .< 0.05) ? true : false
  # Normality test for residuals using the Kolmogorov-Smirnov test
  normal = try
    if n < 100
      ExactOneSampleKSTest(HypothesisTests.ksstats(residual, fit_mle(Normal, residual))...) |> p_result
    else
      ApproximateOneSampleKSTest(HypothesisTests.ksstats(residual, fit_mle(Normal, residual))...) |> p_result
    end
    catch
      false
  end
  # Homoscedasticity test using the White test
  homoscedasticity = try
    WhiteTest([ones(Int(n)) ŷ], residual) |> p_result
  catch
    false
  end
  # Return a matrix containing the calculated parameters
  parameters = [adj_r² syx_in_percentage RMSE MAE AIC coefs_significant normal homoscedasticity]

  return parameters
end

"""
The `criteria_table` function evaluates and ranks multiple regression models based on specified criteria. It generates a comprehensive table of performance metrics for each model, calculates ranks for these metrics, and combines them into a final score. The function allows for flexible selection of evaluation criteria and can return either all models or only the top models based on the combined ranking.

# Parameters:
- `model`: 
    The regression model(s) to be evaluated and compared. This parameter can accept:
    - **Single Linear Regression Model (`TableRegressionModel{<:LinearModel}`)**:
      Evaluates a single linear regression model.
    - **Vector of Linear Regression Models (`Vector{<:TableRegressionModel{<:LinearModel}}`)**:
      Evaluates and compares multiple linear regression models.
    - **Single Grouped linear Regression Model (`GroupedLinearModel`)**:
      Evaluates a regression model that has been grouped by specific factors within a `GroupedLinearModel` object.
    - **Vector of Grouped Linear Regression Models (`Vector{<:GroupedLinearModel}`)**:
      Evaluates and compares multiple grouped regression models.
    
- `criteria::Symbol...`: 
    A variable number of symbols representing the evaluation criteria to include. Possible values include:
    - `:adjr2`: Adjusted R², a measure of the model's explanatory power, adjusted for the number of predictors.
    - `:syx`: Standard error of the estimate (Syx) expressed as a percentage of the mean of the dependent variable (y), indicating the precision of the model's predictions.
    - `:rmse`: Root Mean Squared Error, indicating the average magnitude of residuals.
    - `:mae`: Mean Absolute Error, another accuracy measure based on average absolute residuals.
    - `:aic`: Akaike Information Criterion, balancing goodness of fit with model complexity.
    - `:significance`: Evaluates whether model coefficients are statistically significant.
    - `:normality`: Assesses the normality of residuals using the Kolmogorov-Smirnov test, ensuring that residuals follow a normal distribution.
    - `:homoscedasticity`: Checks for constant variance in residuals using the White test.
    
    If no criteria are specified, the function will use all available criteria by default.

- `best::Union{Bool, Int}=10`: 
    Specifies the number of top models to return based on the combined ranking.
    - `false`: Returns the full table with all models ranked.
    - Integer value: If less than the total number of models, returns only the top `best` models.
  
- `weight::Int=2`: 
    A weighting factor plus `10` applied to the significance, normality, and homoscedasticity criteria to give them more influence in the combined rank.

# Returns:
- `DataFrame`: 
    A sorted table with the evaluated models and their respective metrics. Includes a combined rank based on the selected criteria.

# Examples:
- **Single Model**: 
  `criteria_table(model, :adjr2, :rmse)`
  
- **Vector of Models**:
  `criteria_table([model1, model2], :aic, :mae)`
  
- **Grouped Model**:
  `criteria_table(grouped_model, :significance, :normality)`
  
- **Vector of Grouped Models**:
  `criteria_table([grouped_model1, grouped_model2], :aic, :mae)`
"""
function criteria_table(model::Vector{<:TableRegressionModel{<:LinearModel}}, 
                        criteria::Symbol...; 
                        best::Union{Bool, Int}=10,
                        weight::Int=2) :: DataFrame

  allowed_fields = [:adjr2, :syx, :rmse, :mae, :aic, :significance, :normality, :homoscedasticity]

  if isempty(criteria)
    selected_criteria = allowed_fields
  else
    selected_criteria = Any[s for s in criteria]
  end

  if !issubset(Set(selected_criteria), Set(allowed_fields))
    not_allowed = join(setdiff(selected_criteria, allowed_fields), ", :")
    allowed_msg = "\nAllowed fields are: :" * join(allowed_fields, ", :")
    throw(ArgumentError(":$not_allowed not allowed." * allowed_msg))
  end

  # Generate the criteria parameters for each model in the input vector
  criteria_params = vcat(_criteria_parameters.(model)...)

  # Create a DataFrame from the full criteria parameters
  ct = DataFrame(criteria_params, allowed_fields)

  # Filter the DataFrame columns based on the selected criteria
  ct = select(ct, selected_criteria)

  # Insert the model objects into the DataFrame
  insertcols!(ct, 1, "model" => model)

  # Calculate ranks for each criterion
  ranks = Dict()
  if :adjr2 in selected_criteria
    ranks[:adjr2] = competerank(ct[!, "adjr2"], rev = true)
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
  if :significance in selected_criteria
    ranks[:significance] = competerank(ct[!, "significance"], rev = true) * 10weight
  end
  if :normality in selected_criteria
    ranks[:normality] = competerank(ct[!, "normality"], rev = true) * 10weight
  end
  if :homoscedasticity in selected_criteria
    ranks[:homoscedasticity] = competerank(ct[!, "homoscedasticity"], rev = true) * 10weight
  end

  # Combine ranks into a single score
  combined_rank = sum([ranks[crit] for crit in selected_criteria])

  # Insert the combined rank into the DataFrame
  insertcols!(ct, 2, :rank => combined_rank)

  # Sort the DataFrame by the combined rank
  sort!(ct, :rank)

  # If 'best' is false, return the full DataFrame
  if best === false
    return ct
  elseif best isa Int && best < size(ct, 1)
    # If 'best' is a number and less than the number of models, return the top 'best' models
    top_models = ct[1:best, 1]
    return criteria_table(top_models, criteria...; best = false, weight = weight) # Re-run with selected models
  else
    # Otherwise, return the full DataFrame
    return ct
  end
end

@inline criteria_table(model::TableRegressionModel{<:LinearModel}, criteria::Symbol...; weight::Int=2) :: DataFrame = criteria_table([model], criteria..., weight = weight)

function criteria_table(models::GroupedLinearModel, criteria::Symbol...; weight::Int=2)

  allowed_fields = [:adjr2, :syx, :rmse, :mae, :aic, :significance, :normality, :homoscedasticity]

  if isempty(criteria)
    selected_criteria = allowed_fields
  else
    selected_criteria = Any[s for s in criteria]
  end

  if !issubset(Set(selected_criteria), Set(allowed_fields))
    not_allowed = join(setdiff(selected_criteria, allowed_fields), ", :")
    allowed_msg = "\nAllowed fields are: :" * join(allowed_fields, ", :")
    throw(ArgumentError(":$not_allowed not allowed." * allowed_msg))
  end

  full_model = DataFrame(_criteria_parameters(models), allowed_fields)

  # Filter the DataFrame columns based on the selected criteria
  full_model = select(full_model, selected_criteria)
  
  insertcols!(full_model, 1, "model" => models)
  insertcols!(full_model, 2, :rank => 0)
  insertcols!(full_model, 1, join(models.group_names, " ") => "")

  grouped_models = vcat(models.grouped_models |> values |> collect)
  group_table = DataFrame(:model => grouped_models)
  insertcols!(group_table, 1, join(models.group_names, " ") => collect(keys(models.grouped_models)))
  group_table = innerjoin(group_table, criteria_table(grouped_models, selected_criteria..., best=false, weight=weight), on=:model)


  return vcat(full_model, group_table)
end

function criteria_table(models::Vector{<:GroupedLinearModel{<:LinearModel}}, criteria::Symbol...; weight::Int=2) :: DataFrame
  
  allowed_fields = [:adjr2, :syx, :rmse, :mae, :aic, :significance, :normality, :homoscedasticity]

  if isempty(criteria)
    selected_criteria = allowed_fields
  else
    selected_criteria = Any[s for s in criteria]
  end

  if !issubset(Set(selected_criteria), Set(allowed_fields))
    not_allowed = join(setdiff(selected_criteria, allowed_fields), ", :")
    allowed_msg = "\nAllowed fields are: :" * join(allowed_fields, ", :")
    throw(ArgumentError(":$not_allowed not allowed." * allowed_msg))
  end

  # Generate the criteria parameters for each model in the input vector
  criteria_params = vcat(_criteria_parameters.(models)...)

  # Create a DataFrame from the full criteria parameters
  ct = DataFrame(criteria_params, allowed_fields)

  # Filter the DataFrame columns based on the selected criteria
  ct = select(ct, selected_criteria)

  # Insert the models objects into the DataFrame
  insertcols!(ct, 1, "grouped_models" => models)

  # Calculate ranks for each criterion
  ranks = Dict()
  if :adjr2 in selected_criteria
    ranks[:adjr2] = competerank(ct[!, "adjr2"], rev = true)
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
  if :significance in selected_criteria
    ranks[:significance] = competerank(ct[!, "significance"], rev = true) * 10weight
  end
  if :normality in selected_criteria
    ranks[:normality] = competerank(ct[!, "normality"], rev = true) * 10weight
  end
  if :homoscedasticity in selected_criteria
    ranks[:homoscedasticity] = competerank(ct[!, "homoscedasticity"], rev = true) * 10weight
  end

  # Combine ranks into a single score
  combined_rank = sum([ranks[crit] for crit in selected_criteria])

  # Insert the combined rank into the DataFrame
  insertcols!(ct, 2, :rank => combined_rank)

  # Sort the DataFrame by the combined rank
  sort!(ct, :rank)
  return ct
end

"""
The `criteria_selection` function evaluates and ranks a vector of regression models based on specified criteria, returning the best model according to the combined ranking.

# Parameters:
- `model::Vector{<:TableRegressionModel{<:LinearModel}}`: 
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

- `weight::Int=2`: 
    A weighting factor applied to the significance, normality, and homoscedasticity criteria to give them more influence in the combined rank.

# Returns:
- `TableRegressionModel`: 
    The best model based on the combined ranking of the specified criteria.
"""
@inline criteria_selection(model::Vector{<:TableRegressionModel{<:LinearModel}}, criteria::Symbol...; weight::Int=2) :: TableRegressionModel{<:LinearModel} = criteria_table(model, criteria..., best = 5, weight = weight)[1, 1]