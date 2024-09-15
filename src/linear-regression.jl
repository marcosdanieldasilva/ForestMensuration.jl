function _create_matrix_formulas!(matrix_formulas::Vector{MatrixTerm}, x_term_list::Vector{MixTerm}, q_term::AbstractTerm...)
  for (k, x) ∈ enumerate(x_term_list)
    if isempty(q_term)
      matrix_formulas[k] = MatrixTerm(β0 + x)
    else
      matrix_formulas[k] = MatrixTerm(β0 + x + sum(q_term))
    end
  end
end

function _fit_regression!(fitted_models::Vector{TableRegressionModel}, 
                          y_term_list::Vector{AbstractTerm}, 
                          matrix_formulas::Vector{MatrixTerm}, 
                          cols::NamedTuple,
                          type::Union{Type{LinearModel}, Type{GeneralizedLinearModel}}
                          )

  model_cols = Dict{AbstractTerm, Union{Vector{<:Real}, Nothing}}()
  model_matrix = Dict{MatrixTerm, Union{Matrix{<:Real}, Nothing}}()
  null_schema = Schema()

  for y ∈ y_term_list
    try
      model_cols[y] = modelcols(y, cols)
    catch
      model_cols[y] = nothing
    end
  end

  for x ∈ matrix_formulas
    try
      model_matrix[x] = modelmatrix(x, cols)
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
      # Compute the regression
      if type == LinearModel
        fitted_model = GLM.fit(LinearModel, X, Y)
      else
        fitted_model = GLM.fit(GeneralizedLinearModel, X, Y, Normal(), LogLink())
      end
      # Construct the formula
      formula = FormulaTerm(y, x)
      # Construct the ModelFrame
      mf = ModelFrame(formula, null_schema, cols, type)
      # Construct the ModelMatrix
      mm = ModelMatrix(X, asgn(formula))
      # Pass the fitted_model to TableRegressionModel structure
      push!(fitted_models, TableRegressionModel(fitted_model, mf, mm))
    catch
      # Handle exceptions silently
    end
  end
end

"""
The `regression` function in Julia automatically generates and evaluates multiple simple regression models based on the provided data, including both continuous and categorical variables. This function significantly expands the traditional analysis typically applied in forest biometrics, such as the relationship between tree height and diameter at breast height (DBH), by automatically generating and evaluating 240 unique combinations of dependent and independent variable transformations.

# Parameters:
- `y::Symbol`: 
    The dependent variable (response variable) for the regression model. Typically represents a continuous measure such as tree height or volume in forestry studies.
  
- `x::Symbol`: 
    The independent variable (predictor variable) used to explain variations in the dependent variable. Often represents a measure like diameter at breast height (DBH) in forestry.

- `data::AbstractDataFrame`: 
    The dataset containing the variables for regression. The data frame must include all variables specified in `y`, `x`, and `q`, and it will automatically remove any rows with missing values before performing the regression.

- `q::Symbol...` (optional): 
    A variable number of symbols representing additional categorical variables (qualitative factors) to include in the regression model. These are treated as factors and can influence the model differently based on their levels.

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

- **Standard Regression**:
  The function performs a regression analysis by automatically generating a wide array of possible models. It creates multiple transformations of the dependent and independent variables, combining them into various model forms. The results can be evaluated, and the best models can be selected based on criteria such as adjusted R², RMSE, AIC, and more, using the `criteria_table` function.

- **Qualitative Variables**:
  The function allows the inclusion of categorical variables (`q`) in the regression model. These variables are automatically treated as factors and can be used to capture variations in the dependent variable that are related to these qualitative factors.

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
function regression(y::S, x::S, data::AbstractDataFrame, q::S...; type::Union{Type{LinearModel}, Type{GeneralizedLinearModel}}=LinearModel) where S <: Symbol

  new_data = dropmissing(data[:, [y, x, q...]])

  n, k = size(new_data)

  if n < k + 2
    error("There are not enough data points to perform regression. At least $(k + 2) observations are required.")
  end

  cols = columntable(new_data)

  y_term = concrete_term(term(y), cols, ContinuousTerm)
  x_term = concrete_term(term(x), cols, ContinuousTerm)
  q_term = [concrete_term(term(terms), cols, CategoricalTerm) for terms ∈ q]

  if type == LinearModel
    y_term_list = _dependent_variable(y_term, x_term)
  else
    y_term_list = _dependent_variable(y_term)
  end

  x_term_list = _independent_variable(x_term)

  matrix_formulas = Vector{MatrixTerm}(undef, length(x_term_list))
  _create_matrix_formulas!(matrix_formulas, x_term_list, q_term...)

  fitted_models = Vector{TableRegressionModel}()

  _fit_regression!(fitted_models, y_term_list, matrix_formulas, cols, type)

  isempty(fitted_models) && error("Failed to fit any models")

  return fitted_models
end