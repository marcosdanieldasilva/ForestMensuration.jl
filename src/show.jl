import Base: show

"""
Display the equation of the fitted linear model.

# Arguments
- `io::IO`: The IO stream to which the output is written.
- `fitted_model::FittedLinearModel`: The fitted linear model to be displayed.
"""
function show(io::IO, fitted_model::FittedLinearModel)
  β = coef(fitted_model)
  n = length(β)
  output = string(StatsModels.coefnames(fitted_model.formula.lhs)) * " = $(round(β[1], digits = 4))"
  
  for i in 2:n
    term = coefnames(fitted_model)[i]
    product = string(round(abs(β[i]), sigdigits = 4)) * " * " * term
    output *= signbit(β[i]) ? " - $(product)" : " + $(product)"
  end
  
  print(io, output)
end