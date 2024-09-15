function regression(z::S, y::S, x::S, data::AbstractDataFrame, q::S...) where S <: Symbol

  new_data = dropmissing(data[:, [z, y, x, q...]])

  n, k = size(new_data)

  if n < k + 2
    error("There are not enough data points to perform regression. At least $(k + 2) observations are required.")
  end

  cols = columntable(new_data)

  z_term = concrete_term(term(z), cols, ContinuousTerm)
  y_term = concrete_term(term(y), cols, ContinuousTerm)
  x_term = concrete_term(term(x), cols, ContinuousTerm)
  q_term = [concrete_term(term(terms), cols, CategoricalTerm) for terms ∈ q]

  z_term_list = _dependent_variable(z_term)

  combined_term_list = _generate_combined_terms(y_term, x_term)

  matrix_formulas = Vector{MatrixTerm}(undef, length(combined_term_list))
  _create_matrix_formulas!(matrix_formulas, combined_term_list, q_term...)

  fitted_models = Vector{TableRegressionModel}()

  _fit_regression!(fitted_models, z_term_list, matrix_formulas, cols, GeneralizedLinearModel)

  isempty(fitted_models) && error("Failed to fit any models")

  return fitted_models
end