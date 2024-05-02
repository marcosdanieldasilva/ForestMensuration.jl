import Base: show

function show(io::IO, fitted_model::FittedLinearModel)

  β = coef(fitted_model)
  n = length(β)
  output = StatsModels.coefnames(fitted_model.formula.lhs)  * " = $(round(β[1], digits = 4))"

  for i in 2:n
    terms = coefnames(fitted_model)[i]
    product = string(round(abs(β[i]), sigdigits = 4)) * " * " * terms
    signbit(β[i]) ? output *= " - $(product)" : output *= " + $(product)"
  end

  print(io, output)

end

# function show(io::IO, report::StratifiedReport)
#   println(io)
#   println(io, "ANOVA")
#   show(io, report.anova)
#   println(io)
#   println(io)
#   println(io, "Auxiliary Table")
#   show(io, report.auxiliary_table)
#   println(io)
#   println(io)
#   println(io, "Result Table")
#   show(io, report.result_table)
# end