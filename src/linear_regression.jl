function _fit_linear_model(ft::FormulaTerm, Y::Vector{<:Real}, X::Matrix{<:Real}, data::NamedTuple)
  # Extract the dependent variable (y), independent variables (x), and additional parameters (q) from data
  (y, x, q...) = data
  # Compute the mean of the dependent variable (ȳ) for later calculations
  ȳ = mean(y)
  # Compute the initial product of the transpose of the design matrix and the dependent variable
  β = X'Y
  # Perform in-place Cholesky decomposition on X'X
  chol = cholesky!(X'X)
  # Solve for regression coefficients β in-place using the Cholesky factor
  ldiv!(chol, β)  # β = (X'X)⁻¹ * (X'Y)
  # Allocate space for predicted values (ŷ) with the same structure as Y
  ŷ = similar(Y)
  # Compute predicted values in-place (ŷ = X * β)
  mul!(ŷ, X, β)
  # Determine the number of observations (n) and the number of predictors (ncoef)
  (n, ncoef) = size(X)
  # Calculate the degrees of freedom for residuals
  dof_residuals = n - ncoef
  # Compute the residuals (difference between observed and predicted values)
  residual = Y - ŷ
  # Calculate the Sum of Squared Residuals (SSR)
  SSR = residual ⋅ residual
  # Compute the variance of residuals (σ²) adjusted by degrees of freedom
  σ² = (residual ⋅ residual) / dof_residuals
  # Correct the predicted values and residuals for models with a function on the left-hand side of the formula
  if isa(ft.lhs, FunctionTerm)
    # Adjust predicted values for transformed response
    ŷ = _predict(ft, x, ŷ, σ²)
    # Recompute residuals using adjusted predicted values
    residual = y - ŷ
    SSR = residual ⋅ residual
    # Compute the variance of residuals (σ²) adjusted
    σ²adjusted = SSR / dof_residuals
  else
    σ²adjusted = σ²
  end
  # Calculate the Total Sum of Squares (SST) based on the deviation of y from its mean
  SST = sum(abs2.(y .- ȳ))
  # Compute the coefficient of determination (R²), a measure of model fit
  r² = 1 - SSR / SST
  # Compute the adjusted R², penalized for the number of predictors
  adjr² = 1 - (1 - r²) * (n - 1) / dof_residuals
  # Willmott’s index of agreement (d)
  z .= @. abs(ŷ - ȳ) + abs(y - ȳ)
  d = 1 - SSR / (z ⋅ z)
  # Calculate the Mean Squared Error (MSE) as SSR divided by the number of observations
  MSE = SSR / n
  # Calculate the Root Mean Squared Error (RMSE), a measure of prediction accuracy
  RMSE = √MSE
  # Calculate the Mean Absolute Error (MAE) as the average absolute residual value
  MAE = mean(abs.(residual))
  # Calculate the standard error of the estimate (Syx) as a percentage of the mean response
  Syx = √σ²adjusted / ȳ * 100
  # Compute the log-likelihood of the model for information criteria
  loglike = -n / 2 * (log(2π * σ²adjusted) + 1)
  # Compute the Akaike Information Criterion (AIC) for model quality
  AIC = -2 * loglike + 2 * ncoef
  # Compute the Bayesian Information Criterion (BIC) for penalized model complexity
  BIC = -2 * loglike + log(n) * ncoef
  # Test for normality of residuals using goodness-of-fit
  normality = goodness_of_fit_test(fit_mle(Normal, residual), residual) |> pvalue > 0.05 ? true : false
  # Perform significance testing for coefficients using t-statistics and p-values
  dispersion = rmul!(inv(chol), σ²)
  standard_errors = sqrt.(diag(dispersion))
  t_values = β ./ standard_errors
  p_values = 2 .* ccdf.(TDist(dof_residuals), abs.(t_values))
  # Check if all coefficients are statistically significant at the 0.01 level (99%)
  significance = all(p_values .< 0.01) ? true : false
  # Package the results into a LinearModel structure
  fitted_models = LinearModel(
    ft, data, β, σ², r², adjr², d, MSE, RMSE, MAE, Syx, AIC, BIC, normality, significance
  )
  # Return the fitted model object
  return fitted_models
end

# function regression(ft::FormulaTerm, data::AbstractDataFrame)
#   data = columntable(data)
#   try
#     Y, X = modelcols(ft, data)
#     return _fit_linear_model(ft, Y, X, data)
#   catch
#     ft = apply_schema(ft, StatsModels.schema(data))
#     Y, X = modelcols(ft, data)
#     return _fit_linear_model(ft, Y, X, data)
#   end
# end

"""
    regression(y::Symbol, x::Symbol, data::AbstractDataFrame, q::Symbol...)
  
The `regression` function in Julia automatically generates and evaluates multiple simple regression models based on the provided data, including both continuous and categorical variables. This function significantly expands the traditional analysis typically applied in forest biometrics, such as the relationship between tree height and diameter at breast height (DBH), by automatically generating and evaluating 240 unique combinations of dependent and independent variable transformations.

# Parameters:
- `y::Symbol`: The dependent variable (response variable) for the regression model. Typically represents a continuous measure such as tree height or volume in forestry studies.
  
- `x::Symbol`: The independent variable (predictor variable) used to explain variations in the dependent variable. Often represents a measure like diameter at breast height (DBH) in forestry.

- `data::AbstractDataFrame`: The dataset containing the variables for regression. The data frame must include all variables specified in `y`, `x`, and `q`, and it will automatically remove any rows with missing values before performing the regression.

- `q::Symbol...` (optional): A variable number of symbols representing additional categorical variables (qualitative factors) to include in the regression model. These are treated as factors and can influence the model differently based on their levels.

# Functionality:
- **Dependent Variable Transformations**:
  The function generates multiple transformations of the dependent variable (`y`), including logarithmic and inverse functions. Specifically, the following transformations are created:
  - `log(y)`: Logarithm of `y`.
  - `log_minus(y)`: Logarithm of `y - 1.3`.
  - `log1p(y)`: Logarithm of `1 + y`.
  - `1/y`: Inverse of `y`.
  - `1/(y - 1.3)`: Inverse of `y - 1.3`.
  - `1/√y`: Inverse of the square root of `y`.
  - `1/√(y - 1.3)`: Inverse of the square root of `y - 1.3`.
  - `x/√y`: Ratio of `x` to the square root of `y`.
  - `x/√(y - 1.3)`: Ratio of `x` to the square root of `y - 1.3`.
  - `x²/y`: Square of `x` divided by `y`.
  - `x²/(y - 1.3)`: Square of `x` divided by `y - 1.3`.

- **Independent Variable Transformations**:
  Similarly, multiple transformations of the independent variable (`x`) are created:
  - `x²`: Square of `x`.
  - `log(x)`: Logarithm of `x`.
  - `log(x)²`: Square of the logarithm of `x`.
  - `1/x`: Inverse of `x`.
  - `1/x²`: Inverse of Square of `x`.

- **Combined Model Formulations**:
  The function creates a total of 240 combinations of these transformations by pairing the various forms of `y` and `x`. For example, the function considers combinations such as:
  - `y ~ x + x²`
  - `y ~ x + log(x)`
  - `log(y) ~ log(x) + 1/x`
  - `(x / y) ~ x + log(x)²`
  - And many more.

This comprehensive set of models extends beyond the typical scope of forest biometrics, where usually only a few standard models (around five) are tested. By automatically exploring a wide array of potential relationships, this method allows for a more thorough investigation of the possible functional forms that describe the relationship between tree height and DBH or other relevant variables.

- **Standard Regression**: The function performs a regression analysis by automatically generating a wide array of possible models. It creates multiple transformations of the dependent and independent variables, combining them into various model forms. The results can be evaluated, and the best models can be selected based on criteria such as adjusted R², RMSE, AIC, and more, using the `criteria_table` function.

- **Qualitative Variables**: The function allows the inclusion of categorical variables (`q`) in the regression model. These variables are automatically treated as factors and can be used to capture variations in the dependent variable that are related to these qualitative factors.

# Applications:
This method is ideal for forestry researchers and practitioners who want to ensure they are not overlooking any potentially significant models by expanding their analysis to include a broader range of possible relationships between the variables.

# Examples:
```julia
# Perform standard regression without grouping
models = regression(:height, :diameter, data)

# View the top models based on a specific criteria
best_models = criteria_table(models, :adjr2, :rmse)

```
"""
function regression(y::Symbol, x::Symbol, data::AbstractDataFrame, q::Symbol...)
  # Remover linhas com dados faltantes nas colunas relevantes
  new_data = dropmissing(data[:, [y, x, q...]])

  n, k = size(new_data)

  if n < k + 2
    error("There are not enough data points to perform regression. At least $(k + 2) observations are required.")
  end

  #Attempt to coerce the y and x columns to the Continuous scitype (e.g., Float64)
  try
    d = Dict{Symbol,Type}()
    push!(d, y => ScientificTypes.Continuous, x => ScientificTypes.Continuous)
    map(i -> push!(d, q[i] => ScientificTypes.Multiclass), eachindex(q))
    coerce!(new_data, d)
  catch
    # If coercion fails, an error will be thrown
    error("Unable to coerce variables '$(y)' and '$(x)' to Continuous. Please ensure they contain numeric values.")
  end

  cols = columntable(new_data)

  # Criar termos concretos para Y, X1, X2 e Q
  y_term = concrete_term(term(y), cols, ContinuousTerm)
  x_term = concrete_term(term(x), cols, ContinuousTerm)
  q_terms = [concrete_term(term(qq), cols, CategoricalTerm) for qq in q]

  # Obter a lista de termos dependentes
  y_term_list = _dependent_variable(y_term, x_term)

  # Construir o dicionário para variáveis dependentes
  model_dependent_variable = Dict{AbstractTerm,Union{AbstractVector,Nothing}}()

  for yt in y_term_list
    model_dependent_variable[yt] = try
      modelcols(yt, cols)
    catch
      nothing
    end
  end

  # Preparar armazenamento para as matrizes do modelo
  model_matrix = _independent_variable(x_term, cols, q_terms...)

  fitted_models = Vector{LinearModel}()

  # Solve for each combination of y and matrix formula
  for y in y_term_list
    Y = model_dependent_variable[y]

    if Y === nothing
      continue
    end

    for (mt, X) in model_matrix
      try
        push!(fitted_models, _fit_linear_model(FormulaTerm(y, mt), Y, X, cols))
      catch
      end
    end
  end

  isempty(fitted_models) && error("Failed to fit any models")

  return fitted_models
end

function regression(y::Symbol, x1::Symbol, x2::Symbol, data::AbstractDataFrame, q::Symbol...)
  # Remover linhas com dados faltantes nas colunas relevantes
  new_data = dropmissing(data[:, [y, x1, x2, q...]])

  n, k = size(new_data)

  if n < k + 2
    error("There are not enough data points to perform regression. At least $(k + 2) observations are required.")
  end

  #Attempt to coerce the y and x columns to the Continuous scitype (e.g., Float64)
  try
    d = Dict{Symbol,Type}()
    push!(d, y => ScientificTypes.Continuous, x1 => ScientificTypes.Continuous, x2 => ScientificTypes.Continuous)
    map(i -> push!(d, q[i] => ScientificTypes.Multiclass), eachindex(q))
    coerce!(new_data, d)
  catch
    # If coercion fails, an error will be thrown
    error("Unable to coerce variables '$(y)', '$(x1)' and '$(x2)' to Continuous. Please ensure they contain numeric values.")
  end

  cols = columntable(new_data)

  # Criar termos concretos para Y, X1, X2 e Q
  y_term = concrete_term(term(y), cols, ContinuousTerm)
  x1_term = concrete_term(term(x1), cols, ContinuousTerm)
  x2_term = concrete_term(term(x2), cols, ContinuousTerm)
  q_terms = [concrete_term(term(qq), cols, CategoricalTerm) for qq in q]

  # Obter a lista de termos dependentes
  y_term_list = _dependent_variable(y_term)

  # Construir o dicionário para variáveis dependentes
  model_dependent_variable = Dict{AbstractTerm,Union{AbstractVector,Nothing}}()

  for yt in y_term_list
    model_dependent_variable[yt] = try
      modelcols(yt, cols)
    catch
      nothing
    end
  end

  # Preparar armazenamento para as matrizes do modelo
  model_matrix = _independent_variable(x1_term, x2_term, cols, q_terms...)

  fitted_models = Vector{LinearModel}()

  # Solve for each combination of y and matrix formula
  for y in y_term_list
    Y = model_dependent_variable[y]

    if Y === nothing
      continue
    end

    for (mt, X) in model_matrix
      try
        push!(fitted_models, _fit_linear_model(FormulaTerm(y, mt), Y, X, cols))
      catch
      end
    end
  end

  isempty(fitted_models) && error("Failed to fit any models")

  return fitted_models
end

function grouped_regression(y::Symbol, x::Symbol, data::AbstractDataFrame, g::Symbol...)
  isempty(g) && throw(ArgumentError("At least one grouping variable (g) must be provided after data, e.g., group_regression(:y, :x, data, :group1, :group2)."))
  g = vcat(g...)
  group_data = dropmissing(data[:, [y, x, g...]])
  # Extract group attributes and initialize variables
  groups = RecipesPipeline._extract_group_attributes(tuple(eachcol(group_data[:, g])...))
  labels, idxs = getfield(groups, 1), getfield(groups, 2)
  # Perform general regression
  general_regression = criteria_selection(regression(y, x, data))
  # Perform  qualitative regression
  qualy_regression = criteria_selection(regression(y, x, data, g...), :adjr2, :syx, :aic, :bic, :normality)
  # Perform group-specific regressions
  grouped_models = Dict{String,LinearModel}()

  for (i, label) in enumerate(labels)
    grouped_models[label] = criteria_selection(regression(y, x, group_data[idxs[i], :]))
  end
  # Return the fitted model
  return GroupedLinearModel(general_regression, qualy_regression, grouped_models, g)
end

function grouped_regression(y::Symbol, x1::Symbol, x2::Symbol, data::AbstractDataFrame, g::Symbol...)
  isempty(g) && throw(ArgumentError("At least one grouping variable (g) must be provided."))
  g = vcat(g...)
  group_data = dropmissing(data[:, [y, x1, x2, g...]])
  # Extract group attributes and initialize variables
  groups = RecipesPipeline._extract_group_attributes(tuple(eachcol(group_data[:, g])...))
  labels, idxs = getfield(groups, 1), getfield(groups, 2)
  # Perform general regression
  general_regression = criteria_selection(regression(y, x1, x2, data))
  # Perform  qualitative regression
  qualy_regression = criteria_selection(regression(y, x1, x2, data, g...), :adjr2, :syx, :aic, :bic, :normality)
  # Perform group-specific regressions
  grouped_models = Dict{String,LinearModel}()

  for (i, label) in enumerate(labels)
    grouped_models[label] = criteria_selection(regression(y, x1, x2, group_data[idxs[i], :]))
  end
  # Return the fitted model
  return GroupedLinearModel(general_regression, qualy_regression, grouped_models, g)
end