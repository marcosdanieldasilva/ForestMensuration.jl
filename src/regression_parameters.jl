log_minus(y::Real) = log(y - 1.3)
one_by_y(y::Real) = 1 / y
one_by_y_minus(y::Real) = 1 / (y - 1.3)
one_by_sqrt(y::Real) = 1 / √y
one_by_sqrt_minus(y::Real) = 1 / √(y - 1.3)
x_by_sqrt_y(x::Real, y::Real) = x / √y
x_by_sqrt_y_minus(x::Real, y::Real) = x / √(y - 1.3)
square_x_by_y(x::Real, y::Real) = x^2 / y
square_x_by_y_minus(x::Real, y::Real) = x^2 / (y - 1.3)

# Generates a list of transformed dependent variable terms.
function _dependent_variable(y_term::AbstractTerm, x_term::AbstractTerm)::Vector{AbstractTerm}
  return [
    y_term,
    FunctionTerm(log, [y_term], :(log($(y_term)))),
    FunctionTerm(log_minus, [y_term], :(log_minus($(y_term) - 1.3))),
    FunctionTerm(log1p, [y_term], :(log1p($(y_term)))),
    FunctionTerm(one_by_y, [y_term], :(1 / ($(y_term)))),
    FunctionTerm(one_by_y_minus, [y_term], :(1 / ($(y_term) - 1.3))),
    FunctionTerm(one_by_sqrt, [y_term], :(1 / √($(y_term)))),
    FunctionTerm(one_by_sqrt_minus, [y_term], :(1 / √($(y_term) - 1.3))),
    FunctionTerm(x_by_sqrt_y, [x_term, y_term], :($(x_term) / √$(y_term))),
    FunctionTerm(x_by_sqrt_y_minus, [x_term, y_term], :($(x_term) / √($(y_term) - 1.3))),
    FunctionTerm(square_x_by_y, [x_term, y_term], :($(x_term)^2 / $(y_term))),
    FunctionTerm(square_x_by_y_minus, [x_term, y_term], :($(x_term)^2 / ($(y_term) - 1.3)))
  ]
end

# Generates a list of transformed dependent variable term.
function _dependent_variable(y_term::AbstractTerm)::Vector{AbstractTerm}
  return [
    y_term,
    FunctionTerm(log, [y_term], :(log($(y_term)))),
    FunctionTerm(log1p, [y_term], :(log1p($(y_term))))
  ]
end


# Optimized function to generate combined term matrices
function _independent_variable(x_term::AbstractTerm, cols::NamedTuple, q_terms::AbstractTerm...)
  # Define the six base terms
  x2 = FunctionTerm(x -> x^2, [x_term], :($(x_term)^2))
  log_x = FunctionTerm(log, [x_term], :(log($(x_term))))
  log_x2 = FunctionTerm(x -> log(x)^2, [x_term], :(log($(x_term))^2))
  inv_x = FunctionTerm(x -> 1 / x, [x_term], :($(x_term)^-1))
  inv_x2 = FunctionTerm(x -> 1 / x^2, [x_term], :($(x_term)^-2))
  # Combine all terms into a Dict an calculate de transformed values
  d = Dict{AbstractTerm,Array{Float64}}(t => modelcols(t, cols) for t in [β0, x_term, x2, log_x, log_x2, inv_x, inv_x2])
  # Use the base terms in combinations
  x_terms = MixTerm[
    (x_term,),
    (x2,),
    (log_x,),
    (log_x2,),
    (inv_x,),
    (inv_x2,),
    x_term+x2,
    x_term+log_x,
    x_term+log_x2,
    x_term+inv_x,
    x_term+inv_x2,
    x2+log_x,
    x2+log_x2,
    x2+inv_x,
    x2+inv_x2,
    log_x+log_x2,
    log_x+inv_x,
    log_x+inv_x2,
    log_x2+inv_x,
    log_x2+inv_x2,
    inv_x+inv_x2,
    x_term+x2+log_x,
    x_term+x2+log_x2,
    x_term+x2+inv_x,
    x_term+x2+inv_x2,
    x_term+log_x+log_x2,
    x_term+log_x+inv_x,
    x_term+log_x+inv_x2,
    x_term+log_x2+inv_x,
    x_term+log_x2+inv_x2,
    x_term+inv_x+inv_x2,
    x2+log_x+log_x2,
    x2+log_x+inv_x,
    x2+log_x+inv_x2,
    x2+log_x2+inv_x,
    x2+log_x2+inv_x2,
    x2+inv_x+inv_x2,
    log_x+log_x2+inv_x,
    log_x+log_x2+inv_x2,
    log_x+inv_x+inv_x2,
    log_x2+inv_x+inv_x2
  ]
  # Calculate additional terms if q_terms are provided
  if isempty(q_terms)
    q_sum_term = nothing
  else
    q_sum_term = sum(q_terms)
    q_matrix = modelmatrix(q_terms, cols)
  end
  # Build model_matrix using dictionary comprehension
  model_matrix = Dict{MatrixTerm,Matrix{Float64}}(
    q_sum_term === nothing
    ? MatrixTerm(β0 + x) => hcat(d[β0], [d[t] for t in x]...)
    : MatrixTerm(β0 + x + q_sum_term) => hcat(d[β0], [d[t] for t in x]..., q_matrix)
    for x in x_terms
  )

  return model_matrix

end

function _independent_variable(x1_term::AbstractTerm, x2_term::AbstractTerm, cols::NamedTuple, q_terms::AbstractTerm...)
  # Define transformations for the x1_term variable
  x1_2 = FunctionTerm(x -> x^2, [x1_term], :($(x1_term)^2))
  log_x1 = FunctionTerm(log, [x1_term], :(log($(x1_term))))
  # Define transformations for the x2_term variable 
  x2_2 = FunctionTerm(x -> x^2, [x2_term], :($(x2_term)^2))
  log_x2 = FunctionTerm(log, [x2_term], :(log($(x2_term))))
  # Interactions of two terms
  x1_term_x2_term = x1_term & x2_term
  x1_term_x2_2 = x1_term & x2_2
  x1_term_log_x2 = x1_term & log_x2
  x1_2_x2_term = x1_2 & x2_term
  x1_2_x2_2 = x1_2 & x2_2
  x1_2_log_x2 = x1_2 & log_x2
  log_x1_x2_term = log_x1 & x2_term
  log_x1_x2_2 = log_x1 & x2_2
  log_x1_log_x2 = log_x1 & log_x2
  # Combine all terms into a Dict an calculate de transformed values
  d = Dict{AbstractTerm,Array{Float64}}(
    t => modelcols(t, cols) for t in [
      β0, x1_term, x1_2, log_x1, x2_term, x2_2, log_x2, x1_term_x2_term, x1_term_x2_2,
      x1_term_log_x2, x1_2_x2_term, x1_2_x2_2, x1_2_log_x2, log_x1_x2_term, log_x1_x2_2, log_x1_log_x2
    ]
  )
  # Define x_terms with combinations of two and three terms
  x_terms = MixTerm[
    # Sums of two terms,
    x1_term+x2_term,
    x1_term+x2_2,
    x1_term+log_x2,
    x1_2+x2_term,
    x1_2+x2_2,
    x1_2+log_x2,
    log_x1+x2_term,
    log_x1+x2_2,
    log_x1+log_x2,
    # Interactions of two terms
    (x1_term & x2_term,),
    (x1_term & x2_2,),
    (x1_term & log_x2,),
    (x1_2 & x2_term,),
    (x1_2 & x2_2,),
    (x1_2 & log_x2,),
    (log_x1 & x2_term,),
    (log_x1 & x2_2,),
    (log_x1 & log_x2,),
    x1_term_x2_term+x1_term_x2_2,
    x1_term_x2_term+x1_term_log_x2,
    x1_term_x2_term+x1_2_x2_term,
    x1_term_x2_term+x1_2_x2_2,
    x1_term_x2_term+x1_2_log_x2,
    x1_term_x2_term+log_x1_x2_term,
    x1_term_x2_term+log_x1_x2_2,
    x1_term_x2_term+log_x1_log_x2,
    x1_term_x2_2+x1_term_log_x2,
    x1_term_x2_2+x1_2_x2_term,
    x1_term_x2_2+x1_2_x2_2,
    x1_term_x2_2+x1_2_log_x2,
    x1_term_x2_2+log_x1_x2_term,
    x1_term_x2_2+log_x1_x2_2,
    x1_term_x2_2+log_x1_log_x2,
    x1_term_log_x2+x1_2_x2_term,
    x1_term_log_x2+x1_2_x2_2,
    x1_term_log_x2+x1_2_log_x2,
    x1_term_log_x2+log_x1_x2_term,
    x1_term_log_x2+log_x1_x2_2,
    x1_term_log_x2+log_x1_log_x2,
    x1_2_x2_term+x1_2_x2_2,
    x1_2_x2_term+x1_2_log_x2,
    x1_2_x2_term+log_x1_x2_term,
    x1_2_x2_term+log_x1_x2_2,
    x1_2_x2_term+log_x1_log_x2,
    x1_2_x2_2+x1_2_log_x2,
    x1_2_x2_2+log_x1_x2_term,
    x1_2_x2_2+log_x1_x2_2,
    x1_2_x2_2+log_x1_log_x2,
    x1_2_log_x2+log_x1_x2_term,
    x1_2_log_x2+log_x1_x2_2,
    x1_2_log_x2+log_x1_log_x2,
    log_x1_x2_term+log_x1_x2_2,
    log_x1_x2_term+log_x1_log_x2,
    log_x1_x2_2+log_x1_log_x2,
    # # triple coefs with single interaction 
    x1_term+x2_term+x2_2,
    x1_term+x2_term+log_x2,
    x1_term+x2_term+x1_2,
    x1_term+x2_term+log_x1,
    x1_term+x2_term+x1_term_x2_term,
    x1_term+x2_term+x1_term_x2_2,
    x1_term+x2_term+x1_term_log_x2,
    x1_term+x2_term+x1_2_x2_term,
    x1_term+x2_term+x1_2_x2_2,
    x1_term+x2_term+x1_2_log_x2,
    x1_term+x2_term+log_x1_x2_term,
    x1_term+x2_term+log_x1_x2_2,
    x1_term+x2_term+log_x1_log_x2,
    x1_term+x2_2+log_x2,
    x1_term+x2_2+x1_2,
    x1_term+x2_2+log_x1,
    x1_term+x2_2+x1_term_x2_term,
    x1_term+x2_2+x1_term_x2_2,
    x1_term+x2_2+x1_term_log_x2,
    x1_term+x2_2+x1_2_x2_term,
    x1_term+x2_2+x1_2_x2_2,
    x1_term+x2_2+x1_2_log_x2,
    x1_term+x2_2+log_x1_x2_term,
    x1_term+x2_2+log_x1_x2_2,
    x1_term+x2_2+log_x1_log_x2,
    x1_term+log_x2+x1_2,
    x1_term+log_x2+log_x1,
    x1_term+log_x2+x1_term_x2_term,
    x1_term+log_x2+x1_term_x2_2,
    x1_term+log_x2+x1_term_log_x2,
    x1_term+log_x2+x1_2_x2_term,
    x1_term+log_x2+x1_2_x2_2,
    x1_term+log_x2+x1_2_log_x2,
    x1_term+log_x2+log_x1_x2_term,
    x1_term+log_x2+log_x1_x2_2,
    x1_term+log_x2+log_x1_log_x2,
    x1_2+x2_term+x2_2,
    x1_2+x2_term+log_x2,
    x1_2+x2_term+log_x1,
    x1_2+x2_term+x1_term_x2_term,
    x1_2+x2_term+x1_term_x2_2,
    x1_2+x2_term+x1_term_log_x2,
    x1_2+x2_term+x1_2_x2_term,
    x1_2+x2_term+x1_2_x2_2,
    x1_2+x2_term+x1_2_log_x2,
    x1_2+x2_term+log_x1_x2_term,
    x1_2+x2_term+log_x1_x2_2,
    x1_2+x2_term+log_x1_log_x2,
    x1_2+x2_2+log_x2,
    x1_2+x2_2+log_x1,
    x1_2+x2_2+x1_term_x2_term,
    x1_2+x2_2+x1_term_x2_2,
    x1_2+x2_2+x1_term_log_x2,
    x1_2+x2_2+x1_2_x2_term,
    x1_2+x2_2+x1_2_x2_2,
    x1_2+x2_2+x1_2_log_x2,
    x1_2+x2_2+log_x1_x2_term,
    x1_2+x2_2+log_x1_x2_2,
    x1_2+x2_2+log_x1_log_x2,
    x1_2+log_x2+log_x1,
    x1_2+log_x2+x1_term_x2_term,
    x1_2+log_x2+x1_term_x2_2,
    x1_2+log_x2+x1_term_log_x2,
    x1_2+log_x2+x1_2_x2_term,
    x1_2+log_x2+x1_2_x2_2,
    x1_2+log_x2+x1_2_log_x2,
    x1_2+log_x2+log_x1_x2_term,
    x1_2+log_x2+log_x1_x2_2,
    x1_2+log_x2+log_x1_log_x2,
    log_x1+x2_term+x2_2,
    log_x1+x2_term+log_x2,
    log_x1+x2_term+x1_term_x2_term,
    log_x1+x2_term+x1_term_x2_2,
    log_x1+x2_term+x1_term_log_x2,
    log_x1+x2_term+x1_2_x2_term,
    log_x1+x2_term+x1_2_x2_2,
    log_x1+x2_term+x1_2_log_x2,
    log_x1+x2_term+log_x1_x2_term,
    log_x1+x2_term+log_x1_x2_2,
    log_x1+x2_term+log_x1_log_x2,
    log_x1+x2_2+log_x2,
    log_x1+x2_2+x1_term_x2_term,
    log_x1+x2_2+x1_term_x2_2,
    log_x1+x2_2+x1_term_log_x2,
    log_x1+x2_2+x1_2_x2_term,
    log_x1+x2_2+x1_2_x2_2,
    log_x1+x2_2+x1_2_log_x2,
    log_x1+x2_2+log_x1_x2_term,
    log_x1+x2_2+log_x1_x2_2,
    log_x1+x2_2+log_x1_log_x2,
    log_x1+log_x2+x1_term_x2_term,
    log_x1+log_x2+x1_term_x2_2,
    log_x1+log_x2+x1_term_log_x2,
    log_x1+log_x2+x1_2_x2_term,
    log_x1+log_x2+x1_2_x2_2,
    log_x1+log_x2+x1_2_log_x2,
    log_x1+log_x2+log_x1_x2_term,
    log_x1+log_x2+log_x1_x2_2,
    log_x1+log_x2+log_x1_log_x2,
    # triple coefs with doble interaction 
    x1_term+x1_term_x2_term+x1_term_x2_2,
    x1_term+x1_term_x2_term+x1_term_log_x2,
    x1_term+x1_term_x2_term+x1_2_x2_term,
    x1_term+x1_term_x2_term+x1_2_x2_2,
    x1_term+x1_term_x2_term+x1_2_log_x2,
    x1_term+x1_term_x2_term+log_x1_x2_term,
    x1_term+x1_term_x2_term+log_x1_x2_2,
    x1_term+x1_term_x2_term+log_x1_log_x2,
    x1_term+x1_term_x2_2+x1_term_log_x2,
    x1_term+x1_term_x2_2+x1_2_x2_term,
    x1_term+x1_term_x2_2+x1_2_x2_2,
    x1_term+x1_term_x2_2+x1_2_log_x2,
    x1_term+x1_term_x2_2+log_x1_x2_term,
    x1_term+x1_term_x2_2+log_x1_x2_2,
    x1_term+x1_term_x2_2+log_x1_log_x2,
    x1_term+x1_term_log_x2+x1_2_x2_term,
    x1_term+x1_term_log_x2+x1_2_x2_2,
    x1_term+x1_term_log_x2+x1_2_log_x2,
    x1_term+x1_term_log_x2+log_x1_x2_term,
    x1_term+x1_term_log_x2+log_x1_x2_2,
    x1_term+x1_term_log_x2+log_x1_log_x2,
    x1_term+x1_2_x2_term+x1_2_x2_2,
    x1_term+x1_2_x2_term+x1_2_log_x2,
    x1_term+x1_2_x2_term+log_x1_x2_term,
    x1_term+x1_2_x2_term+log_x1_x2_2,
    x1_term+x1_2_x2_term+log_x1_log_x2,
    x1_term+x1_2_x2_2+x1_2_log_x2,
    x1_term+x1_2_x2_2+log_x1_x2_term,
    x1_term+x1_2_x2_2+log_x1_x2_2,
    x1_term+x1_2_x2_2+log_x1_log_x2,
    x1_term+x1_2_log_x2+log_x1_x2_term,
    x1_term+x1_2_log_x2+log_x1_x2_2,
    x1_term+x1_2_log_x2+log_x1_log_x2,
    x1_term+log_x1_x2_term+log_x1_x2_2,
    x1_term+log_x1_x2_term+log_x1_log_x2,
    x1_term+log_x1_x2_2+log_x1_log_x2,
    x1_2+x1_term_x2_term+x1_term_x2_2,
    x1_2+x1_term_x2_term+x1_term_log_x2,
    x1_2+x1_term_x2_term+x1_2_x2_term,
    x1_2+x1_term_x2_term+x1_2_x2_2,
    x1_2+x1_term_x2_term+x1_2_log_x2,
    x1_2+x1_term_x2_term+log_x1_x2_term,
    x1_2+x1_term_x2_term+log_x1_x2_2,
    x1_2+x1_term_x2_term+log_x1_log_x2,
    x1_2+x1_term_x2_2+x1_term_log_x2,
    x1_2+x1_term_x2_2+x1_2_x2_term,
    x1_2+x1_term_x2_2+x1_2_x2_2,
    x1_2+x1_term_x2_2+x1_2_log_x2,
    x1_2+x1_term_x2_2+log_x1_x2_term,
    x1_2+x1_term_x2_2+log_x1_x2_2,
    x1_2+x1_term_x2_2+log_x1_log_x2,
    x1_2+x1_term_log_x2+x1_2_x2_term,
    x1_2+x1_term_log_x2+x1_2_x2_2,
    x1_2+x1_term_log_x2+x1_2_log_x2,
    x1_2+x1_term_log_x2+log_x1_x2_term,
    x1_2+x1_term_log_x2+log_x1_x2_2,
    x1_2+x1_term_log_x2+log_x1_log_x2,
    x1_2+x1_2_x2_term+x1_2_x2_2,
    x1_2+x1_2_x2_term+x1_2_log_x2,
    x1_2+x1_2_x2_term+log_x1_x2_term,
    x1_2+x1_2_x2_term+log_x1_x2_2,
    x1_2+x1_2_x2_term+log_x1_log_x2,
    x1_2+x1_2_x2_2+x1_2_log_x2,
    x1_2+x1_2_x2_2+log_x1_x2_term,
    x1_2+x1_2_x2_2+log_x1_x2_2,
    x1_2+x1_2_x2_2+log_x1_log_x2,
    x1_2+x1_2_log_x2+log_x1_x2_term,
    x1_2+x1_2_log_x2+log_x1_x2_2,
    x1_2+x1_2_log_x2+log_x1_log_x2,
    x1_2+log_x1_x2_term+log_x1_x2_2,
    x1_2+log_x1_x2_term+log_x1_log_x2,
    x1_2+log_x1_x2_2+log_x1_log_x2,
    x2_term+x1_term_x2_term+x1_term_x2_2,
    x2_term+x1_term_x2_term+x1_term_log_x2,
    x2_term+x1_term_x2_term+x1_2_x2_term,
    x2_term+x1_term_x2_term+x1_2_x2_2,
    x2_term+x1_term_x2_term+x1_2_log_x2,
    x2_term+x1_term_x2_term+log_x1_x2_term,
    x2_term+x1_term_x2_term+log_x1_x2_2,
    x2_term+x1_term_x2_term+log_x1_log_x2,
    x2_term+x1_term_x2_2+x1_term_log_x2,
    x2_term+x1_term_x2_2+x1_2_x2_term,
    x2_term+x1_term_x2_2+x1_2_x2_2,
    x2_term+x1_term_x2_2+x1_2_log_x2,
    x2_term+x1_term_x2_2+log_x1_x2_term,
    x2_term+x1_term_x2_2+log_x1_x2_2,
    x2_term+x1_term_x2_2+log_x1_log_x2,
    x2_term+x1_term_log_x2+x1_2_x2_term,
    x2_term+x1_term_log_x2+x1_2_x2_2,
    x2_term+x1_term_log_x2+x1_2_log_x2,
    x2_term+x1_term_log_x2+log_x1_x2_term,
    x2_term+x1_term_log_x2+log_x1_x2_2,
    x2_term+x1_term_log_x2+log_x1_log_x2,
    x2_term+x1_2_x2_term+x1_2_x2_2,
    x2_term+x1_2_x2_term+x1_2_log_x2,
    x2_term+x1_2_x2_term+log_x1_x2_term,
    x2_term+x1_2_x2_term+log_x1_x2_2,
    x2_term+x1_2_x2_term+log_x1_log_x2,
    x2_term+x1_2_x2_2+x1_2_log_x2,
    x2_term+x1_2_x2_2+log_x1_x2_term,
    x2_term+x1_2_x2_2+log_x1_x2_2,
    x2_term+x1_2_x2_2+log_x1_log_x2,
    x2_term+x1_2_log_x2+log_x1_x2_term,
    x2_term+x1_2_log_x2+log_x1_x2_2,
    x2_term+x1_2_log_x2+log_x1_log_x2,
    x2_term+log_x1_x2_term+log_x1_x2_2,
    x2_term+log_x1_x2_term+log_x1_log_x2,
    x2_term+log_x1_x2_2+log_x1_log_x2,
    x2_2+x1_term_x2_term+x1_term_x2_2,
    x2_2+x1_term_x2_term+x1_term_log_x2,
    x2_2+x1_term_x2_term+x1_2_x2_term,
    x2_2+x1_term_x2_term+x1_2_x2_2,
    x2_2+x1_term_x2_term+x1_2_log_x2,
    x2_2+x1_term_x2_term+log_x1_x2_term,
    x2_2+x1_term_x2_term+log_x1_x2_2,
    x2_2+x1_term_x2_term+log_x1_log_x2,
    x2_2+x1_term_x2_2+x1_term_log_x2,
    x2_2+x1_term_x2_2+x1_2_x2_term,
    x2_2+x1_term_x2_2+x1_2_x2_2,
    x2_2+x1_term_x2_2+x1_2_log_x2,
    x2_2+x1_term_x2_2+log_x1_x2_term,
    x2_2+x1_term_x2_2+log_x1_x2_2,
    x2_2+x1_term_x2_2+log_x1_log_x2,
    x2_2+x1_term_log_x2+x1_2_x2_term,
    x2_2+x1_term_log_x2+x1_2_x2_2,
    x2_2+x1_term_log_x2+x1_2_log_x2,
    x2_2+x1_term_log_x2+log_x1_x2_term,
    x2_2+x1_term_log_x2+log_x1_x2_2,
    x2_2+x1_term_log_x2+log_x1_log_x2,
    x2_2+x1_2_x2_term+x1_2_x2_2,
    x2_2+x1_2_x2_term+x1_2_log_x2,
    x2_2+x1_2_x2_term+log_x1_x2_term,
    x2_2+x1_2_x2_term+log_x1_x2_2,
    x2_2+x1_2_x2_term+log_x1_log_x2,
    x2_2+x1_2_x2_2+x1_2_log_x2,
    x2_2+x1_2_x2_2+log_x1_x2_term,
    x2_2+x1_2_x2_2+log_x1_x2_2,
    x2_2+x1_2_x2_2+log_x1_log_x2,
    x2_2+x1_2_log_x2+log_x1_x2_term,
    x2_2+x1_2_log_x2+log_x1_x2_2,
    x2_2+x1_2_log_x2+log_x1_log_x2,
    x2_2+log_x1_x2_term+log_x1_x2_2,
    x2_2+log_x1_x2_term+log_x1_log_x2,
    x2_2+log_x1_x2_2+log_x1_log_x2,
    log_x2+x1_term_x2_term+x1_term_x2_2,
    log_x2+x1_term_x2_term+x1_term_log_x2,
    log_x2+x1_term_x2_term+x1_2_x2_term,
    log_x2+x1_term_x2_term+x1_2_x2_2,
    log_x2+x1_term_x2_term+x1_2_log_x2,
    log_x2+x1_term_x2_term+log_x1_x2_term,
    log_x2+x1_term_x2_term+log_x1_x2_2,
    log_x2+x1_term_x2_term+log_x1_log_x2,
    log_x2+x1_term_x2_2+x1_term_log_x2,
    log_x2+x1_term_x2_2+x1_2_x2_term,
    log_x2+x1_term_x2_2+x1_2_x2_2,
    log_x2+x1_term_x2_2+x1_2_log_x2,
    log_x2+x1_term_x2_2+log_x1_x2_term,
    log_x2+x1_term_x2_2+log_x1_x2_2,
    log_x2+x1_term_x2_2+log_x1_log_x2,
    log_x2+x1_term_log_x2+x1_2_x2_term,
    log_x2+x1_term_log_x2+x1_2_x2_2,
    log_x2+x1_term_log_x2+x1_2_log_x2,
    log_x2+x1_term_log_x2+log_x1_x2_term,
    log_x2+x1_term_log_x2+log_x1_x2_2,
    log_x2+x1_term_log_x2+log_x1_log_x2,
    log_x2+x1_2_x2_term+x1_2_x2_2,
    log_x2+x1_2_x2_term+x1_2_log_x2,
    log_x2+x1_2_x2_term+log_x1_x2_term,
    log_x2+x1_2_x2_term+log_x1_x2_2,
    log_x2+x1_2_x2_term+log_x1_log_x2,
    log_x2+x1_2_x2_2+x1_2_log_x2,
    log_x2+x1_2_x2_2+log_x1_x2_term,
    log_x2+x1_2_x2_2+log_x1_x2_2,
    log_x2+x1_2_x2_2+log_x1_log_x2,
    log_x2+x1_2_log_x2+log_x1_x2_term,
    log_x2+x1_2_log_x2+log_x1_x2_2,
    log_x2+x1_2_log_x2+log_x1_log_x2,
    log_x2+log_x1_x2_term+log_x1_x2_2,
    log_x2+log_x1_x2_term+log_x1_log_x2,
    log_x2+log_x1_x2_2+log_x1_log_x2,
  ]
  # Calculate additional terms if q_terms are provided
  if isempty(q_terms)
    q_sum_term = nothing
  else
    q_sum_term = sum(q_terms)
    q_matrix = modelmatrix(q_terms, cols)
  end
  # Build model_matrix using dictionary comprehension
  model_matrix = Dict{MatrixTerm,Matrix{Float64}}(
    q_sum_term === nothing
    ? MatrixTerm(β0 + x) => hcat(d[β0], [d[t] for t in x]...)
    : MatrixTerm(β0 + x + q_sum_term) => hcat(d[β0], [d[t] for t in x]..., q_matrix)
    for x in x_terms
  )

  return model_matrix
end