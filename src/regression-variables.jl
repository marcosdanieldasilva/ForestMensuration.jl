@inline log_minus(y::Real) :: Float64 = log(y - 1.3)
@inline log_inverse(y::Real) :: Float64 = log(1 / y)
@inline log_inverse_minus(y::Real) :: Float64 = log(1 / (y - 1.3))
@inline one_by_y(y::Real) :: Float64 = 1 / y
@inline one_by_y_minus(y::Real) :: Float64 = 1 / (y - 1.3)
@inline one_by_sqrt(y::Real) :: Float64 = 1 / √y
@inline one_by_sqrt_minus(y::Real) :: Float64 = 1 / √(y - 1.3)
@inline one_by_cbrt(y::Real) :: Float64 = 1 / ∛y
@inline one_by_cbrt_minus(y::Real) :: Float64 = 1 / ∛(y - 1.3)
@inline x_by_y(x::Real, y::Real) :: Float64 = x / y
@inline x_by_y_minus(x::Real, y::Real) :: Float64 = x / (y - 1.3)
@inline x_by_sqrt_y(x::Real, y::Real) :: Float64 = x / √y
@inline x_by_sqrt_y_minus(x::Real, y::Real) :: Float64 = x / √(y - 1.3)
@inline x_by_cbrt_y(x::Real, y::Real) :: Float64 = x / ∛y
@inline x_by_cbrt_y_minus(x::Real, y::Real) :: Float64 = x / ∛(y - 1.3)
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
    FunctionTerm(log_minus, [y_term], :(log($(y_term) - 1.3))),
    FunctionTerm(log_inverse, [y_term], :(log($(y_term) ^- 1))),
    FunctionTerm(log_inverse_minus, [y_term], :(log(($(y_term) - 1.3) ^- 1))),
    FunctionTerm(one_by_y, [y_term], :(1/($(y_term)))),
    FunctionTerm(one_by_y_minus, [y_term], :(1/($(y_term) - 1.3))),
    FunctionTerm(one_by_sqrt, [y_term], :(1/√($(y_term)))),
    FunctionTerm(one_by_sqrt_minus, [y_term], :(1/√($(y_term) - 1.3))),
    FunctionTerm(one_by_cbrt, [y_term], :(1/∛($(y_term)))),
    FunctionTerm(one_by_cbrt_minus, [y_term], :(1/∛($(y_term) - 1.3))),
    # FunctionTerm(x_by_y, [x_term, y_term], :($(x_term) / $(y_term)))
    # FunctionTerm(x_by_y_minus, [x_term, y_term], :($(x_term) / ($(y_term) - 1.3)))
    FunctionTerm(x_by_sqrt_y, [x_term, y_term], :($(x_term) / √$(y_term))),
    FunctionTerm(x_by_sqrt_y_minus, [x_term, y_term], :($(x_term) / √($(y_term) - 1.3))),
    # FunctionTerm(x_by_cbrt_y, [x_term, y_term], :($(x_term) / ∛$(y_term)))
    # FunctionTerm(x_by_cbrt_y_minus, [x_term, y_term], :($(x_term) / ∛($(y_term) - 1.3)))
    FunctionTerm(square_x_by_y, [x_term, y_term], :($(x_term) ^ 2 / $(y_term))),
    FunctionTerm(square_x_by_y_minus, [x_term, y_term], :($(x_term) ^ 2 / ($(y_term) - 1.3)))
  ]
end

function _dependent_variable(y_term::AbstractTerm) :: Vector{AbstractTerm}
  return [
    y_term,
    FunctionTerm(log, [y_term], :(log($(y_term)))),
    FunctionTerm(log_inverse, [y_term], :(log($(y_term) ^- 1)))
    # FunctionTerm(one_by_y, [y_term], :(1/($(y_term)))),
    # FunctionTerm(one_by_sqrt, [y_term], :(1/√($(y_term)))),
    # FunctionTerm(one_by_cbrt, [y_term], :(1/∛($(y_term))))
    ]
end

"""
Generates a list of transformed independent variable terms and their interactions.

# Arguments
- `x_term::AbstractTerm`: The independent variable term.

# Returns
- `Vector{MixTerm}`: A list of transformed independent variable terms and their interactions.
"""
# function _independent_variable(x_term::AbstractTerm) :: Vector{MixTerm}
#   x_terms = [
#     x_term
#     FunctionTerm(x -> x ^ 2, [x_term], :($(x_term) ^ 2))
#     FunctionTerm(log, [x_term], :(log($(x_term))))
#     FunctionTerm(x -> log(x) ^ 2, [x_term], :(log($(x_term)) ^ 2))
#     FunctionTerm(x -> log(1 / x), [x_term], :(log(1 / $(x_term))))
#     FunctionTerm(x -> 1 / x, [x_term], :($(x_term) ^- 1))
#     FunctionTerm(x -> 1 / x ^ 2, [x_term], :($(x_term) ^- 2))
#     FunctionTerm(x -> log(x / (1 + x)), [x_term], :(log($(x_term) / (1 + $(x_term)))))
#   ]

#   n = length(x_terms)
#   x_term_list = Vector{MixTerm}(undef, binomial(n, 2))
#   k = 1

#   @inbounds for i in 1:n-1
#     for j in i+1:n
#       x_term_list[k] = x_terms[i] + x_terms[j]
#       k += 1
#     end
#   end

#   return x_term_list
# end
function _independent_variable(x_term::AbstractTerm) :: Vector{MixTerm}
  return [
    x_term,
    FunctionTerm(x -> x ^ 2, [x_term], :($(x_term) ^ 2)),
    FunctionTerm(log, [x_term], :(log($(x_term)))),
    FunctionTerm(x -> log(x) ^ 2, [x_term], :(log($(x_term)) ^ 2)),
    # FunctionTerm(x -> log(1 / x), [x_term], :(log(1 / $(x_term)))),
    FunctionTerm(x -> 1 / x, [x_term], :($(x_term) ^ -1)),
    FunctionTerm(x -> 1 / x ^ 2, [x_term], :($(x_term) ^ -2)),

    x_term + FunctionTerm(x -> x ^ 2, [x_term], :($(x_term) ^ 2)),
    x_term + FunctionTerm(log, [x_term], :(log($(x_term)))),
    x_term + FunctionTerm(x -> log(x) ^ 2, [x_term], :(log($(x_term)) ^ 2)),
    # x_term + FunctionTerm(x -> log(1 / x), [x_term], :(log(1 / $(x_term)))),
    x_term + FunctionTerm(x -> 1 / x, [x_term], :($(x_term) ^ -1)),
    x_term + FunctionTerm(x -> 1 / x ^ 2, [x_term], :($(x_term) ^ -2)),

    FunctionTerm(x -> x ^ 2, [x_term], :($(x_term) ^ 2)) + FunctionTerm(log, [x_term], :(log($(x_term)))),
    FunctionTerm(x -> x ^ 2, [x_term], :($(x_term) ^ 2)) + FunctionTerm(x -> log(x) ^ 2, [x_term], :(log($(x_term)) ^ 2)),
    # FunctionTerm(x -> x ^ 2, [x_term], :($(x_term) ^ 2)) + FunctionTerm(x -> log(1 / x), [x_term], :(log(1 / $(x_term)))),
    FunctionTerm(x -> x ^ 2, [x_term], :($(x_term) ^ 2)) + FunctionTerm(x -> 1 / x, [x_term], :($(x_term) ^ -1)),

    FunctionTerm(log, [x_term], :(log($(x_term)))) + FunctionTerm(x -> log(x) ^ 2, [x_term], :(log($(x_term)) ^ 2)),
    # FunctionTerm(log, [x_term], :(log($(x_term)))) + FunctionTerm(x -> log(1 / x), [x_term], :(log(1 / $(x_term)))),
    FunctionTerm(log, [x_term], :(log($(x_term)))) + FunctionTerm(x -> 1 / x, [x_term], :($(x_term) ^ -1)),
    FunctionTerm(log, [x_term], :(log($(x_term)))) + FunctionTerm(x -> 1 / x ^ 2, [x_term], :($(x_term) ^ -2)),

    # FunctionTerm(x -> log(x) ^ 2, [x_term], :(log($(x_term)) ^ 2)) + FunctionTerm(x -> log(1 / x), [x_term], :(log(1 / $(x_term)))),
    FunctionTerm(x -> log(x) ^ 2, [x_term], :(log($(x_term)) ^ 2)) + FunctionTerm(x -> 1 / x, [x_term], :($(x_term) ^ -1)),
    FunctionTerm(x -> log(x) ^ 2, [x_term], :(log($(x_term)) ^ 2)) + FunctionTerm(x -> 1 / x ^ 2, [x_term], :($(x_term) ^ -2)),

    # FunctionTerm(x -> log(1 / x), [x_term], :(log(1 / $(x_term)))) + FunctionTerm(x -> 1 / x, [x_term], :($(x_term) ^ -1)),
    # FunctionTerm(x -> log(1 / x), [x_term], :(log(1 / $(x_term)))) + FunctionTerm(x -> 1 / x ^ 2, [x_term], :($(x_term) ^ -2)),

    FunctionTerm(x -> 1 / x, [x_term], :($(x_term) ^ -1)) + FunctionTerm(x -> 1 / x ^ 2, [x_term], :($(x_term) ^ -2))
  ]
end

# function _independent_variable(x1_term::AbstractTerm, x2_term::AbstractTerm) :: Vector{MixTerm}
#   x1_terms = _independent_variable(x1_term)
#   x2_terms = _independent_variable(x2_term)

#   n = length(x1_terms)^2

#   combined_terms = Vector{MixTerm}(undef, n)

#   k = 1

#   @inbounds for i in 1:n-1
#     for j in i+1:n
#       combined_terms[k] = x1_terms[i] + x2_terms[j]
#       k += 1
#     end
#   end

#   return combined_terms
# end
function _independent_variable(x1_term::AbstractTerm, x2_term::AbstractTerm) :: Vector{MixTerm}
  x1_terms = _independent_variable(x1_term)
  x2_terms = _independent_variable(x2_term)

  n1 = length(x1_terms)
  n2 = length(x2_terms)
  n = n1 * n2

  combined_terms = Vector{MixTerm}(undef, n)

  k = 1
  for i in 1:n1
    for j in 1:n2
      combined_terms[k] = x1_terms[i] + x2_terms[j]
      k += 1
    end
  end

  return combined_terms
end

function generate_combined_terms(x1_term::AbstractTerm, x2_term::AbstractTerm) :: Vector{MixTerm}

  x1_terms = _independent_variable(x1_term)
  x2_terms = _independent_variable(x2_term)

  combined_terms = Vector{MixTerm}()

  for x1 in x1_terms[1:6]
    for x2 in x2_terms
      push!(combined_terms, x1 + x2)
      push!(combined_terms, x1 * x2)
    end
  end

  for x2 in x2_terms[1:6]
    for x1 in x1_terms
      push!(combined_terms, x1 + x2)
      push!(combined_terms, x1 * x2)
    end
  end

  return combined_terms
end