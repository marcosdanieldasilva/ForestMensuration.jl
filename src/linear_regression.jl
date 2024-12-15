function _fit_regression!(
  fitted_models::Vector{TableRegressionModel},
  y_term_list::Vector{AbstractTerm},
  model_matrix::Dict{MatrixTerm,Matrix{Float64}},
  cols::NamedTuple)

  model_dependent_variable = Dict{AbstractTerm,Union{AbstractVector,Nothing}}()

  for yt in y_term_list
    model_dependent_variable[yt] = try
      modelcols(yt, cols)
    catch
      nothing
    end
  end

  for y in y_term_list
    Y = model_dependent_variable[y]

    if Y === nothing
      continue
    end

    for (mt, X) in model_matrix
      try
        # Perform linear regression using GLM.fit
        fitted_model = GLM.fit(LinearModel, X, Y)

        # Create a formula combining y and x terms
        formula = FormulaTerm(y, mt)

        # Create a ModelFrame object to store schema and columns for the model
        mf = ModelFrame(formula, emptySchema, cols, LinearModel)

        # Create a ModelMatrix object based on the formula
        mm = ModelMatrix(X, asgn(formula))

        # Store the fitted model along with its associated data in TableRegressionModel
        push!(fitted_models, TableRegressionModel(fitted_model, mf, mm))
      catch
        # Handle any errors silently during model fitting
      end
    end
  end
end
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
  # remove empty values from the data
  new_data = dropmissing(data[:, [y, x, q...]])

  n, k = size(new_data)

  if n < k + 2
    error("There are not enough data points to perform regression. At least $(k + 2) observations are required.")
  end

  #Attempt to coerce the y and x columns to the Continuous scitype (e.g., Float64)
  try
    coerce!(new_data, y => ScientificTypes.Continuous, x => ScientificTypes.Continuous)
  catch
    # If coercion fails, an error will be thrown
    error("Unable to coerce variables '$(y)' and '$(x)' to Continuous. Please ensure they contain numeric values.")
  end

  # Convert the DataFrame to a column table (named tuple of vectors)
  cols = columntable(new_data)

  # create terms to performe the linear regression
  y_term = concrete_term(term(y), cols, ContinuousTerm)
  x_term = concrete_term(term(x), cols, ContinuousTerm)
  q_terms = [concrete_term(term(qq), cols, CategoricalTerm) for qq in q]

  y_term_list = _dependent_variable(y_term, x_term)

  model_matrix = _independent_variable(x_term, cols, q_terms...)

  fitted_models = Vector{TableRegressionModel}()

  _fit_regression!(fitted_models, y_term_list, model_matrix, cols)

  isempty(fitted_models) && error("Failed to fit any models")

  return fitted_models
end

function regression(y::Symbol, x1::Symbol, x2::Symbol, data::AbstractDataFrame, q::Symbol...)
  # remove empty values from the data
  new_data = dropmissing(data[:, [y, x1, x2, q...]])

  n, k = size(new_data)

  if n < k + 2
    error("There are not enough data points to perform regression. At least $(k + 2) observations are required.")
  end

  #Attempt to coerce the y and x columns to the Continuous scitype (e.g., Float64)
  try
    coerce!(new_data, y => ScientificTypes.Continuous, x1 => ScientificTypes.Continuous, x2 => ScientificTypes.Continuous)
  catch
    # If coercion fails, an error will be thrown
    error("Unable to coerce variables '$(y)', '$(x1)' and '$(x2)' to Continuous. Please ensure they contain numeric values.")
  end

  # Convert the DataFrame to a column table (named tuple of vectors)
  cols = columntable(new_data)

  # create terms to performe the linear regression
  y_term = concrete_term(term(y), cols, ContinuousTerm)
  x1_term = concrete_term(term(x1), cols, ContinuousTerm)
  x2_term = concrete_term(term(x2), cols, ContinuousTerm)
  q_terms = [concrete_term(term(qq), cols, CategoricalTerm) for qq in q]

  y_term_list = _dependent_variable(y_term)

  model_matrix = _independent_variable(x1_term, x2_term, cols, q_terms...)

  fitted_models = Vector{TableRegressionModel}()

  _fit_regression!(fitted_models, y_term_list, model_matrix, cols)

  isempty(fitted_models) && error("Failed to fit any models")

  return fitted_models
end