function _predict(ŷ::Vector{<:Real}, x::Vector{<:Real}, σ²::Real, function_name::Symbol)

  if function_name == :log
    meyer_factor = exp(σ² * 0.5)
    return @. exp(ŷ) * meyer_factor

  elseif function_name == :log_minus
    meyer_factor = exp(σ² * 0.5)
    return @. exp(ŷ) * meyer_factor + 1.3

  elseif function_name == :log_inverse
    meyer_factor = exp(σ² * 0.5)
    return @. 1 / (exp(ŷ) * meyer_factor)

  elseif function_name == :log_inverse_minus
    meyer_factor = exp(σ² * 0.5)
    return @. 1 / (exp(ŷ) * meyer_factor) + 1.3
    
  elseif function_name == :one_by_y
    return @. 1 / ŷ
    
  elseif function_name == :one_by_y_minus
    return @. 1 / ŷ + 1.3

  elseif function_name == :one_by_sqrt
    return @. 1 / ŷ ^ 2

  elseif function_name == :one_by_sqrt_minus
    return @. 1 / ŷ ^ 2 + 1.3

  elseif function_name == :one_by_cbrt
    return @. 1 / ŷ ^ 3

  elseif function_name == :one_by_cbrt_minus
    return @. 1 / ŷ ^ 3 + 1.3

  elseif function_name == :x_by_y
    return @. x / ŷ

  elseif function_name == :x_by_y_minus
    return @. x / ŷ + 1.3

  elseif function_name == :x_by_sqrt_y
    return @. (x / ŷ) ^ 2

  elseif function_name == :x_by_sqrt_y_minus
    return @. (x / ŷ) ^ 2 + 1.3

  elseif function_name == :x_by_cbrt_y
    return @. (x / ŷ) ^ 3

  elseif function_name == :x_by_cbrt_y_minus
    return @. (x / ŷ) ^ 3 + 1.3

  elseif function_name == :square_x_by_y
    return @. x ^ 2 / ŷ

  elseif function_name == :square_x_by_y_minus
    return @. x ^ 2 / ŷ + 1.3

  else
    error("The $function_name don't exist")
  end

end 

function _create_matrix_formulas!(matrix_formulas::Vector{MatrixTerm}, x_term_list::Vector{MixTerm}, q_term::AbstractTerm...; effect::Symbol)
  if isempty(q_term)
    for (k, x) ∈ enumerate(x_term_list)
      matrix_formulas[k] = MatrixTerm(β0 + x)
    end
  else
    n = length(q_term)
    for (k, x) ∈ enumerate(x_term_list)
      if effect == :additive
        sum_term = x
        for i ∈ 1:n
          sum_term += q_term[i]
        end
        matrix_formulas[k] = MatrixTerm(β0 + sum_term)
      else
        prod_term = x
        for i ∈ 1:n
          prod_term *= q_term[i]
        end
        matrix_formulas[k] = MatrixTerm(β0 + prod_term)
      end
    end
  end
end

function _fit_regression!(fitted_models::Vector{FittedLinearModel}, y_term_list::Vector{AbstractTerm}, 
                          matrix_formulas::Vector{MatrixTerm}, data::AbstractDataFrame,
                          y_observed::Vector{<:Real}, x_observed::Vector{<:Real})
  n = size(data, 1)
  ŷ = similar(Vector{Float64}, n)
  model_cols = Dict{AbstractTerm, Union{Vector{<:Real}, Nothing}}()
  model_matrix = Dict{MatrixTerm, Union{Matrix{<:Real}, Nothing}}()

  for y ∈ y_term_list
    try
      push!(model_cols, y => modelcols(y, data))
    catch
      push!(model_cols, y => nothing)
    end
  end

  for x ∈ matrix_formulas
    try
      push!(model_matrix, x => modelmatrix(x, data))
    catch
      push!(model_matrix, x => nothing)
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
    end
  end

end

function regression(y::S, x::S, data::AbstractDataFrame, q::S...; best::Union{Bool, Int}=true, effect::S=:additive) where S <: Symbol

  if best < 0
    throw(ArgumentError("best must be a positive integer or false"))
  elseif  !isempty(q)
    effect ∈ (:additive, :interactive) || throw(ArgumentError("Invalid effect argument. Expected :additive or :interactive."))
  end

  new_data = dropmissing(data[:, [y, x, q...]])

  n, k = size(new_data)

  if n < k + 2
    error("There are not enough data points to perform regression. At least $(k + 2) observations are required.")
  end

  (y_observed, x_observed) = eachcol(new_data)
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