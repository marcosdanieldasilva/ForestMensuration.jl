function _predict(ft::FormulaTerm, x::Vector{<:Real}, ŷ::Vector{<:Real}, σ²::Real)
  function_name = nameof(ft.lhs.f)
  # The Meyer factor is derived from the model's residual variance (σ²) and is used for adjustment
  meyer_factor = exp(σ² / 2)
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

function _predict(model::LinearModel, x::Vector{<:Real}, ŷ::Vector{<:Real})
  function_name = nameof(model.formula.lhs.f)
  # The Meyer factor is derived from the model's residual variance (σ²) and is used for adjustment
  meyer_factor = exp(model.σ² / 2)
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
    predict(model::LinearModel)
  
The `predict` function family provides a versatile way to generate predictions from regression models, 
  supporting both individual and grouped models. It handles predictions on the original scale even if the
   dependent variable (`y`) has been transformed (e.g., `log(y)`), ensuring that any transformations 
   applied during model fitting are correctly reversed, including the application of Meyer correction 
   factors for logarithmic transformations.

# Parameters:
- `model`: 
    The regression model(s) to be evaluated and compared. This parameter can accept:
    - `LinearModel`: A single linear regression model.

# Returns:
- `Vector{<:Real}` or `Vector{Union{Missing, <:Real}}`: The predicted values on the original scale of `y`, adjusted for any transformations and corrected using the Meyer factor for logarithmic transformations.

# Key Features:
- **Handles Transformed Dependent Variables:** If the dependent variable was transformed (e.g., using log transformations), the function correctly inverts the transformation to return predictions on the original scale.
- **Applies Meyer Correction Factor:** For models using logarithmic transformations, the Meyer correction factor is applied to the predictions to correct for the bias introduced by the log transformation.

# Examples:
- **Single Model Prediction:**
  ```julia
  y_pred = predict(model)
  ```
"""
function predict(model::LinearModel)
  # Extract response vector (Y) and predictor matrix (X) from the model
  Y, X = modelcols(model.formula, model.data)
  # Number of observations
  n = length(Y)
  # Allocate memory for predicted values
  ŷ = similar(Vector{Float64}, n)
  # Compute predicted values: ŷ = X * β
  mul!(ŷ, X, model.β)  # In-place multiplication for efficiency
  # Handle special cases where the left-hand side (lhs) is a function term
  if isa(model.formula.lhs, FunctionTerm)
    # Apply the function-specific prediction logic
    ŷ = _predict(model, model.data[2], ŷ)
  end
  # Return the vector of predicted values
  return ŷ
end

"""
    predict(model::LinearModel, data)
Predicts the response variable for a given dataset based on the provided regression model.

# Arguments:

- `model`: 
    The regression model(s) to be evaluated and compared. This parameter can accept:
    - `LinearModel`: A single linear regression model.
    - `GroupedLinearModel`: A grouped linear model where different regression models are fitted for different subsets of data.
- `data`: A dataset compatible with the Tables.jl interface. Must include the predictors required by the model.

# Returns:
- Predicted values as a vector. If there are missing values in the predictors, a predict array is returned with missing values handled.

# Examples:
- **Single Model Prediction:**
  ```julia
  y_pred = predict(model, data)
  ```
"""
function predict(model::LinearModel, data)
  # Ensure the provided data is a valid table or DataFrameRow
  if !(Tables.istable(data) || data isa DataFrameRow)
    throw(ArgumentError("The provided data must be a valid table. Ensure it is a supported type (e.g., DataFrame, DataFrameRow, NamedTuple, or any Tables.jl compatible structure)."))
  end
  # Convert AbstractDataFrame to a column table for compatibility with StatsModels
  if data isa AbstractDataFrame
    data = columntable(data)
  end
  # Convert DataFrameRow to a NamedTuple for compatibility
  if data isa DataFrameRow
    data = columntable(DataFrame(data))
  end
  # Remove missing values and prepare the input data for the model
  x, nonmissings = StatsModels.missing_omit(data, model.formula.rhs)
  # Generate the model matrix (design matrix) from the input data
  X = modelmatrix(model.formula.rhs, x)
  # Compute the predicted values: ŷ = X * β
  ŷ = X * model.β
  # Handle special cases when the formula's left-hand side (lhs) is a function term
  if isa(model.formula.lhs, FunctionTerm)
    ŷ = _predict(model, x[1], ŷ)
  end
  # Return predictions, handling missing values if necessary
  return length(ŷ) == 1 ? ŷ[1] : StatsModels._return_predictions(
    Tables.materializer(data),
    ŷ,
    nonmissings,
    length(nonmissings)
  )
end

"""
    predict!(model::LinearModel, data::AbstractDataFrame)

The `predict!` function computes predictions from a regression model and adds these predictions 
  directly to the provided data frame as new columns. It is particularly useful in forest inventory 
  data where not all trees have been measured for a specific variable, allowing the model to estimate 
    these missing values.

# Parameters:
- `model`: 
    The regression model(s) to be evaluated and compared. This parameter can accept:
    - `LinearModel`: A single linear regression model.

- `data`: The data frame (`AbstractDataFrame`) containing the input data. The function will add new columns to this data frame.

# Functionality:
- **Predicted Values Column (`_predict`)**:
  The function calculates the predicted values for the dependent variable (`y`) based on the input model
     and appends these values as a new column in the data frame with the suffix `_predict`.

- **Real or Estimated Values Column (`_real`)**:
  The function also creates a `_real` column where the actual measured values of `y` are preserved if 
    they exist. For observations where `y` is missing (or set to `0.0`), the predicted value from the 
    model is used instead.

  This setup is ideal for forest inventory datasets where certain tree attributes (like height or volume)
     might not be measured for every tree, and predictions need to be filled in for these gaps.

# Examples:
```julia
# Apply predictions to a data frame
predict!(model, forest_inventory_data)
```
"""
function predict!(model::LinearModel, data::AbstractDataFrame)
  y = propertynames(model.data)[1]
  col_names = Symbol(string(y, "_predict")), Symbol(string(y, "_real"))
  insertcols!(data, col_names[1] => predict(model, data), makeunique=true)
  insertcols!(data, col_names[2] => coalesce.(replace(data[!, y], 0.0 => missing), data[!, col_names[1]]), makeunique=true)
end

"""
    residuals(model::LinearModel)
"""
residuals(model::LinearModel) = model.residuals