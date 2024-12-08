function predict(model::FittedLinearModel)
  Y, X = modelcols(model.formula, model.data)
  # Number of observations and predictor variables
  (n, p) = size(X)
  # Degrees of freedom for residuals
  dof_residuals = n - p
  ŷ = similar(Vector{Float64}, n)
  # Calculate the original predicted values
  mul!(ŷ, X, model.β)

  if isa(model.formula.lhs, FunctionTerm)
    function_name = nameof(model.formula.lhs.f)
    #The Meyer factor is a metric derived from the sum of squared residuals and the degrees of freedom 
    #of the model, which assesses the adjusted variance of the residuals.
    if function_name in [:log, :log_minus, :log1p]
      # Calculate the residuals values
      residual = Y - ŷ
      # Variance of residuals: σ² = Σ(residual^2) / degrees of freedom for residuals
      σ² = (residual ⋅ residual) / dof_residuals
      meyer_factor = exp(σ² / 2)
      ŷ = _predict(function_name, model.data[2], ŷ, meyer_factor=meyer_factor)
    else
      ŷ = _predict(function_name, model.data[2], ŷ)
    end
  end

  return ŷ
end

function _predict(function_name::Symbol, x::Vector{<:Real}, ŷ::Vector{<:Real}; meyer_factor::Real=1.0)
  return @. begin
    if function_name == :log
      exp(ŷ) * meyer_factor
    elseif function_name == :log_minus
      exp(ŷ) * meyer_factor + 1.3
    elseif function_name == :log1p
      expm1(ŷ) * meyer_factor
    elseif function_name == :one_by_y
      1 / ŷ
    elseif function_name == :one_by_y_minus
      1 / ŷ + 1.3
    elseif function_name == :one_by_sqrt
      1 / ŷ^2
    elseif function_name == :one_by_sqrt_minus
      1 / ŷ^2 + 1.3
    elseif function_name == :x_by_sqrt_y
      (x / ŷ)^2
    elseif function_name == :x_by_sqrt_y_minus
      (x / ŷ)^2 + 1.3
    elseif function_name == :square_x_by_y
      x^2 / ŷ
    elseif function_name == :square_x_by_y_minus
      x^2 / ŷ + 1.3
    else
      ŷ
    end
  end
end


"""
Computes the Meyer factor for a fitted linear model.

The Meyer factor is a metric derived from the sum of squared residuals and the degrees of freedom 
of the model, which assesses the adjusted variance of the residuals. It can be useful for model evaluation.

# Arguments:
- `model::FittedLinearModel`: A fitted linear model.

# Returns:
- `meyer_factor::Float64`: The computed Meyer factor.

"""
function meyer_factor(model::FittedLinearModel)
  Y, X = modelcols(f, model.data)
  # Number of observations and predictor variables
  (n, p) = size(X)
  # Degrees of freedom for residuals
  dof_residuals = n - p
  ŷ = similar(Vector{Float64}, n)
  # Calculate the original predicted values
  mul!(ŷ, X, model.β)
  # Calculate the residuals values
  residual = Y - ŷ
  σ² = (residual ⋅ residual) / dof_residuals
  meyer_factor = exp(σ² / 2)

  return meyer_factor
end

function predict(model::FittedLinearModel, data)
  if !Tables.istable(data)
    error()
  end
  if data isa AbstractDataFrame
    data = columntable(data)
  end
  x, nonmissings = StatsModels.missing_omit(data, model.formula.rhs)
  X = modelmatrix(model.formula.rhs, x)
  ŷ = X * model.β

  if isa(model.formula.lhs, FunctionTerm)
    function_name = nameof(model.formula.lhs.f)
    if function_name in [:log, :log_minus, :log1p]
      ŷ = _predict(function_name, model.data[2], ŷ, meyer_factor=meyer_factor(model))
    else
      ŷ = _predict(function_name, model.data[2], ŷ)
    end
  end

  length(unique(nonmissings)) == 1 ? ŷ : StatsModels._return_predictions(Tables.materializer(data), ŷ, nonmissings, length(nonmissings))
end



function _prediction(model::FittedLinearModel, x::Vector, ŷ::Vector)

  f = formula(model)
  function_name = nameof(f.lhs.f)

  if function_name in [:log, :log_minus, :log1p]
    y = model.model.rr.y
    resid = y - predict(model)
    σ² = (resid ⋅ resid) / dof_residual(model)
    meyer_factor = exp(σ² / 2)
  end

  return @. begin
    if function_name == :log
      exp(ŷ) * meyer_factor
    elseif function_name == :log_minus
      exp(ŷ) * meyer_factor + 1.3
    elseif function_name == :log1p
      expm1(ŷ) * meyer_factor
    elseif function_name == :one_by_y
      1 / ŷ
    elseif function_name == :one_by_y_minus
      1 / ŷ + 1.3
    elseif function_name == :one_by_sqrt
      1 / ŷ^2
    elseif function_name == :one_by_sqrt_minus
      1 / ŷ^2 + 1.3
    elseif function_name == :x_by_sqrt_y
      (x / ŷ)^2
    elseif function_name == :x_by_sqrt_y_minus
      (x / ŷ)^2 + 1.3
    elseif function_name == :square_x_by_y
      x^2 / ŷ
    elseif function_name == :square_x_by_y_minus
      x^2 / ŷ + 1.3
    else
      ŷ
    end
  end
end

function _prediction(model::TableRegressionModel{<:GeneralizedLinearModel}, ŷ::Vector)

  resid = model.model.rr.wrkresid
  σ² = (resid ⋅ resid) / dof_residual(model)
  meyer_factor = exp(σ² / 2)

  return ŷ * meyer_factor

end

"""
    prediction(model::TableRegressionModel)
    prediction(model::TableRegressionModel, data::AbstractDataFrame)
    prediction(model::TableRegressionModel, data::DataFrameRow)
  
The `prediction` function family provides a versatile way to generate predictions from regression models, 
  supporting both individual and grouped models. It handles predictions on the original scale even if the
   dependent variable (`y`) has been transformed (e.g., `log(y)`), ensuring that any transformations 
   applied during model fitting are correctly reversed, including the application of Meyer correction 
   factors for logarithmic transformations.

# Parameters:
- `model`: A single linear regression model.

- `data`: 
  The input data for which predictions are needed. This parameter can be:
  - `AbstractDataFrame`: A data frame containing the input data.
  - `DataFrameRow`: A single row of data from a data frame.

# Returns:
- `Vector{<:Real}` or `Vector{Union{Missing, <:Real}}`: The predicted values on the original scale of `y`, adjusted for any transformations and corrected using the Meyer factor for logarithmic transformations.

# Key Features:
- **Handles Transformed Dependent Variables:** If the dependent variable was transformed (e.g., using log transformations), the function correctly inverts the transformation to return predictions on the original scale.
- **Applies Meyer Correction Factor:** For models using logarithmic transformations, the Meyer correction factor is applied to the predictions to correct for the bias introduced by the log transformation.
- **Supports Grouped Models:** When using grouped models, the function automatically selects and applies the correct model for each data point based on its group membership.

# Examples:
- **Single Model Prediction:**
  ```julia
  y_pred = prediction(model, data)
  ```
"""
prediction(model::TableRegressionModel) = return model.model isa GeneralizedLinearModel ? _prediction(model, predict(model)) : model.mf.f.lhs isa FunctionTerm ? _prediction(model, model.mf.data[2], predict(model)) : predict(model)

function prediction(model::TableRegressionModel, data::AbstractDataFrame)
  if model.model isa GeneralizedLinearModel
    return _prediction(model, predict(model, data))
  elseif model.mf.f.lhs isa FunctionTerm
    x, nonmissings = missing_omit(columntable(data), model.mf.f.rhs)
    return _prediction(model, x[1], predict(model, data))
  else
    return predict(model, data)
  end
end

prediction(model::TableRegressionModel, data::DataFrameRow) = prediction(model, DataFrame(data))

"""
    prediction!(model::TableRegressionModel, data::AbstractDataFrame)

The `prediction!` function computes predictions from a regression model and adds these predictions 
  directly to the provided data frame as new columns. It is particularly useful in forest inventory 
  data where not all trees have been measured for a specific variable, allowing the model to estimate 
    these missing values.

# Parameters:
- `model`: A single linear regression model.

- `data`: The data frame (`AbstractDataFrame`) containing the input data. The function will add new columns to this data frame.

# Functionality:
- **Predicted Values Column (`_prediction`)**:
  The function calculates the predicted values for the dependent variable (`y`) based on the input model
     and appends these values as a new column in the data frame with the suffix `_prediction`.

- **Real or Estimated Values Column (`_real`)**:
  The function also creates a `_real` column where the actual measured values of `y` are preserved if 
    they exist. For observations where `y` is missing (or set to `0.0`), the predicted value from the 
    model is used instead.

  This setup is ideal for forest inventory datasets where certain tree attributes (like height or volume)
     might not be measured for every tree, and predictions need to be filled in for these gaps.

# Examples:
```julia
# Apply predictions to a data frame
prediction!(model, forest_inventory_data)
```
"""
function prediction!(model::TableRegressionModel, data::AbstractDataFrame)
  y = propertynames(model.mf.data)[1]
  col_names = Symbol(string(y, "_prediction")), Symbol(string(y, "_real"))
  insertcols!(data, col_names[1] => prediction(model, data), makeunique=true)
  insertcols!(data, col_names[2] => coalesce.(replace(data[!, y], 0.0 => missing), data[!, col_names[1]]), makeunique=true)
end