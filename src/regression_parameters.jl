@inline log_minus(y::Real) :: Float64 = log(y - 1.3)

@inline one_by_y(y::Real) :: Float64 = 1 / y
@inline one_by_y_minus(y::Real) :: Float64 = 1 / (y - 1.3)
@inline one_by_sqrt(y::Real) :: Float64 = 1 / √y
@inline one_by_sqrt_minus(y::Real) :: Float64 = 1 / √(y - 1.3)

@inline x_by_sqrt_y(x::Real, y::Real) :: Float64 = x / √y
@inline x_by_sqrt_y_minus(x::Real, y::Real) :: Float64 = x / √(y - 1.3)
@inline square_x_by_y(x::Real, y::Real) :: Float64 = x ^ 2 / y
@inline square_x_by_y_minus(x::Real, y::Real) :: Float64 = x ^ 2 / (y - 1.3)

"""
Generates a list of transformed dependent variable terms.

# Arguments
- `y_term::AbstractTerm`: The dependent variable term.
- `x_term::AbstractTerm`: The independent variable term.

# Returns
- `Vector{AbstractTerm}`: A list of transformed dependent variable terms.
"""
function _dependent_variable(y_term::AbstractTerm, x_term::AbstractTerm) :: Vector{AbstractTerm}
  return [
    y_term,

    FunctionTerm(log, [y_term], :(log($(y_term)))),
    FunctionTerm(log_minus, [y_term], :(log_minus($(y_term) - 1.3))),
    FunctionTerm(log1p, [y_term], :(log1p($(y_term)))),

    FunctionTerm(one_by_y, [y_term], :(1/($(y_term)))),
    FunctionTerm(one_by_y_minus, [y_term], :(1/($(y_term) - 1.3))),
    FunctionTerm(one_by_sqrt, [y_term], :(1/√($(y_term)))),
    FunctionTerm(one_by_sqrt_minus, [y_term], :(1/√($(y_term) - 1.3))),

    FunctionTerm(x_by_sqrt_y, [x_term, y_term], :($(x_term) / √$(y_term))),
    FunctionTerm(x_by_sqrt_y_minus, [x_term, y_term], :($(x_term) / √($(y_term) - 1.3))),
    FunctionTerm(square_x_by_y, [x_term, y_term], :($(x_term) ^ 2 / $(y_term))),
    FunctionTerm(square_x_by_y_minus, [x_term, y_term], :($(x_term) ^ 2 / ($(y_term) - 1.3)))
  ]
end

"""
Generates a list of transformed dependent variable term.

# Arguments
- `y_term::AbstractTerm`: The dependent variable term.

# Returns
- `Vector{AbstractTerm}`: A list of transformed dependent variable term.
"""
function _dependent_variable(y_term::AbstractTerm) :: Vector{AbstractTerm}
  return [
    y_term,

    FunctionTerm(log, [y_term], :(log($(y_term)))),
    FunctionTerm(log_minus, [y_term], :(log_minus($(y_term) - 1.3))),
    FunctionTerm(log1p, [y_term], :(log1p($(y_term)))),

    FunctionTerm(one_by_y, [y_term], :(1/($(y_term)))),
    FunctionTerm(one_by_y_minus, [y_term], :(1/($(y_term) - 1.3))),
    FunctionTerm(one_by_sqrt, [y_term], :(1/√($(y_term)))),
    FunctionTerm(one_by_sqrt_minus, [y_term], :(1/√($(y_term) - 1.3)))
  ]
end

"""
Generates a list of transformed independent variable terms and their interactions.

# Arguments
- `x_term::AbstractTerm`: The independent variable term.

# Returns
- `Vector{MixTerm}`: A list of transformed independent variable terms and their interactions.
"""
function _independent_variable(x_term::AbstractTerm) :: Vector{MixTerm}
  # Define the six base terms
  x2 = FunctionTerm(x -> x ^ 2, [x_term], :($(x_term) ^ 2))
  log_x = FunctionTerm(log, [x_term], :(log($(x_term))))
  log_x2 = FunctionTerm(x -> log(x) ^ 2, [x_term], :(log($(x_term)) ^ 2))
  inv_x = FunctionTerm(x -> 1 / x, [x_term], :($(x_term) ^ -1))
  inv_x2 = FunctionTerm(x -> 1 / x^2, [x_term], :($(x_term) ^ -2))
  # Use the base terms in combinations
  return [
      x_term,
      x2,
      log_x,
      log_x2,
      inv_x,
      inv_x2,

      x_term + x2,
      x_term + log_x,
      x_term + log_x2,
      x_term + inv_x,
      x_term + inv_x2,

      x2 + log_x,
      x2 + log_x2,
      x2 + inv_x,

      log_x + log_x2,
      log_x + inv_x,
      log_x + inv_x2,

      log_x2 + inv_x,
      log_x2 + inv_x2,

      inv_x + inv_x2
  ]
end

function _generate_combined_terms(x_term::AbstractTerm, y_term::AbstractTerm) :: Vector{MixTerm}
  
  # Define transformations for the x_term variable
  # x^2: square of the x_term
  # log1p: log(1 + x_term) to prevent log(0) issues
  x2 = FunctionTerm(x -> x ^ 2, [x_term], :($(x_term) ^ 2))
  log_x = FunctionTerm(log1p, [x_term], :(log1p($(x_term))))

  # Collect all x transformations into a list
  x_terms = [
    x_term, 
    x2,
    log_x
  ]

  # Define transformations for the y_term variable 
  # y^2: square of the y_term
  # log1p: log(1 + y_term) to prevent log(0) issues
  y2 = FunctionTerm(y -> y ^ 2, [y_term], :($(y_term) ^ 2))
  log_y = FunctionTerm(log1p, [y_term], :(log1p($(y_term))))

  # Collect all y transformations into a list
  y_terms = [
    y_term,
    y2,
    log_y
  ]

  # Generate all possible sums between x_terms and y_terms
  sum_terms = [x + y for x in x_terms for y in y_terms]

  # Generate all possible interactions (products) between x_terms and y_terms
  interaction_terms = [x & y for x in x_terms for y in y_terms]

  # Combine the sum_terms and interaction_terms into a new set of terms
  combined_terms = [st + it for st in sum_terms for it in interaction_terms]

  # The number of interaction terms generated
  n = length(interaction_terms)

  # Loop over x_terms and combine them with pairs of interaction terms
  for x in x_terms
    for i in 1:n-1
      for j in i+1:n
        push!(combined_terms, x + interaction_terms[i] + interaction_terms[j])
      end
    end
  end

  # Similarly, loop over y_terms and combine them with pairs of interaction terms
  for y in y_terms
    for i in 1:n-1
      for j in i+1:n
        push!(combined_terms, y + interaction_terms[i] + interaction_terms[j])
      end
    end
  end

  # Return the final vector of all combined terms
  return combined_terms
end