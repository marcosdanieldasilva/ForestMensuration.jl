"""
Predict values based on different transformation functions.

# Arguments
- `ŷ::Vector{<:Real}`: Predicted values.
- `x::Vector{<:Real}`: Independent variable values.
- `σ²::Real`: Variance of residuals.
- `function_name::Symbol`: Name of the transformation function.

# Returns
- `Vector{<:Real}`: Transformed predicted values.
"""
function _predict(ŷ::Vector{<:Real}, x::Vector{<:Real}, σ²::Real, function_name::Symbol)
  meyer_factor = exp(σ² * 0.5)
  
  return @. begin
    if function_name == :log
      exp(ŷ) * meyer_factor
    elseif function_name == :log_minus
      exp(ŷ) * meyer_factor + 1.3
    elseif function_name == :log_inverse
      1 / (exp(ŷ) * meyer_factor)
    elseif function_name == :log_inverse_minus
      1 / (exp(ŷ) * meyer_factor) + 1.3
    elseif function_name == :one_by_y
      1 / ŷ
    elseif function_name == :one_by_y_minus
      1 / ŷ + 1.3
    elseif function_name == :one_by_sqrt
      1 / ŷ^2
    elseif function_name == :one_by_sqrt_minus
      1 / ŷ^2 + 1.3
    elseif function_name == :one_by_cbrt
      1 / ŷ^3
    elseif function_name == :one_by_cbrt_minus
      1 / ŷ^3 + 1.3
    elseif function_name == :x_by_y
      x / ŷ
    elseif function_name == :x_by_y_minus
      x / ŷ + 1.3
    elseif function_name == :x_by_sqrt_y
      (x / ŷ)^2
    elseif function_name == :x_by_sqrt_y_minus
      (x / ŷ)^2 + 1.3
    elseif function_name == :x_by_cbrt_y
      (x / ŷ)^3
    elseif function_name == :x_by_cbrt_y_minus
      (x / ŷ)^3 + 1.3
    elseif function_name == :square_x_by_y
      x^2 / ŷ
    elseif function_name == :square_x_by_y_minus
      x^2 / ŷ + 1.3
    else
      error("The $function_name don't exist")
    end
  end
end


"""
Create matrix formulas for linear regression models.

# Arguments
- `matrix_formulas::Vector{MatrixTerm}`: Vector to store the resulting matrix formulas.
- `x_term_list::Vector{MixTerm}`: List of independent variable terms.
- `q_term::AbstractTerm...`: Additional terms for the model.
- `effect::Symbol`: Type of effect to consider (`:additive` or `:interactive`).

# Returns
- `Void`: Modifies `matrix_formulas` in place.
"""
function _create_matrix_formulas!(matrix_formulas::Vector{MatrixTerm}, x_term_list::Vector{MixTerm}, q_term::AbstractTerm...; effect::Symbol)
  n = length(q_term)
  for (k, x) ∈ enumerate(x_term_list)
  term = x
  if n > 0
    term = if effect == :additive
    x + sum(q_term)
    else
    x * prod(q_term)
    end
  end
  matrix_formulas[k] = MatrixTerm(β0 + term)
  end
end

"""
Fit linear regression models to the data.

# Arguments
- `fitted_models::Vector{FittedLinearModel}`: Vector to store the fitted models.
- `y_term_list::Vector{AbstractTerm}`: List of dependent variable terms.
- `matrix_formulas::Vector{MatrixTerm}`: List of matrix formulas for the models.
- `data::AbstractDataFrame`: The data frame containing the data.
- `y_observed::Vector{<:Real}`: Observed values of the dependent variable.
- `x_observed::Vector{<:Real}`: Observed values of the independent variable.

# Returns
- `Void`: Modifies `fitted_models` in place.
"""
function _fit_regression!(fitted_models::Vector{FittedLinearModel}, y_term_list::Vector{AbstractTerm}, 
        matrix_formulas::Vector{MatrixTerm}, data::AbstractDataFrame,
        y_observed::Vector{<:Real}, x_observed::Vector{<:Real})
  n = size(data, 1)
  ŷ = similar(Vector{Float64}, n)
  model_cols = Dict{AbstractTerm, Union{Vector{<:Real}, Nothing}}()
  model_matrix = Dict{MatrixTerm, Union{Matrix{<:Real}, Nothing}}()

  for y ∈ y_term_list
  try
    model_cols[y] = modelcols(y, data)
  catch
    model_cols[y] = nothing
  end
  end

  for x ∈ matrix_formulas
  try
    model_matrix[x] = modelmatrix(x, data)
  catch
    model_matrix[x] = nothing
  end
  end

  for (y, x) ∈ Iterators.product(y_term_list, matrix_formulas)
  Y = model_cols[y]
  X = model_matrix[x]

  if Y === nothing || X === nothing
    continue
  end

  try
    # Calculate the regression coefficients
    β = X'Y 
    # Compute the Cholesky decomposition of X'X for optimization
    chol = cholesky!(X'X) 
    # Calculate the coefficients of the fitted regression
    ldiv!(chol, β)
    # Calculate the predicted values
    mul!(ŷ, X, β)
    # Calculate the residuals values
    residual = Y - ŷ
    # Number of predictor variables
    p = size(X, 2)
    # Degrees of freedom for residuals
    dof_residuals = n - p
    # Estimate the variance of residuals
    σ² = (residual ⋅ residual) / dof_residuals 
    # Calculate residuals according to the regression FunctionTerm, if applicable
    if isa(y, FunctionTerm)
    residual = y_observed - _predict(ŷ, x_observed, σ², nameof(y.f))
    end
    # Calculate the RMSE (Root Mean Squared Error)
    RMSE = √((residual ⋅ residual) / n)
    # fit the FormulaTerm
    formula = FormulaTerm(y, x)
    # Pass the fitted model to FittedLinearModel structure
    push!(fitted_models, FittedLinearModel(formula, data, β, σ², RMSE, chol))
  catch
    # Handle exceptions silently
  end
  end
end

"""
Perform regression analysis.

# Arguments
- `y::S`: Symbol representing the dependent variable.
- `x::S`: Symbol representing the independent variable.
- `data::AbstractDataFrame`: Data frame containing the data.
- `q::S...`: Symbols representing additional terms for the model.
- `best::Union{Bool, Int}`: Whether to return the best model(s) based on RMSE.
- `effect::S`: Type of effect to consider (`:additive` or `:interactive`).

# Returns
- `Vector{FittedLinearModel}`: Vector of fitted models.
"""
function regression(y::S, x::S, data::AbstractDataFrame, q::S...; best::Union{Bool, Int}=true, effect::S=:additive) where S <: Symbol
  if best < 0
  throw(ArgumentError("best must be a positive integer or false"))
  elseif !isempty(q)
  effect ∈ (:additive, :interactive) || throw(ArgumentError("Invalid effect argument. Expected :additive or :interactive."))
  end

  new_data = dropmissing(data[:, [y, x, q...]])

  n, k = size(new_data)
  if n < k + 2
  error("There are not enough data points to perform regression. At least $(k + 2) observations are required.")
  end

  y_observed = new_data[!, y]
  x_observed = new_data[!, x]
  y_term = concrete_term(term(y), y_observed, ContinuousTerm)
  x_term = concrete_term(term(x), x_observed, ContinuousTerm)
  q_term = [concrete_term(term(terms), new_data[!, terms], CategoricalTerm) for terms ∈ q]
  y_term_list = _dependent_variable(y_term, x_term)
  x_term_list = _indepedent_variable(x_term)

  matrix_formulas = Vector{MatrixTerm}(undef, length(x_term_list))
  _create_matrix_formulas!(matrix_formulas, x_term_list, q_term..., effect = effect)

  fitted_models = Vector{FittedLinearModel}()
  _fit_regression!(fitted_models, y_term_list, matrix_formulas, new_data, y_observed, x_observed)

  isempty(fitted_models) && error("Failed to fit any models")

  if best == 0
  return fitted_models
  else
  return reduce(vcat, partialsort!(fitted_models, 1:best, by = x -> x.RMSE))
  end
end