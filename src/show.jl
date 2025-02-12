import Base: show

# Custom show method for SiteAnalysis to display  the site_table and site_plot
function show(io::IO, analysis::SiteAnalysis)
  show(io, analysis.site_table)
  if !haskey(ENV, "CI")
    display(analysis.site_plot)
  end
end

_coefnames(model::LinearModel) = vcat("β0", string.(StatsModels.coefnames(model.formula.rhs.terms[2:end])))

# Display the equation of the fitted linear model.
function show(io::IO, model::LinearModel)
  β = model.β
  n = length(β)
  output = string(StatsModels.coefnames(model.formula.lhs)) * " = $(round(β[1], digits = 4))"

  for i in 2:n
    term = _coefnames(model)[i]
    product = string(round(abs(β[i]), sigdigits=4)) * " * " * term
    output *= signbit(β[i]) ? " - $(product)" : " + $(product)"
  end

  print(io, output)
end

function show(io::IO, test::T) where {T<:KolmogorovSmirnovTest}
  println(io, name_of_test(test))
  println(io, repeat("-", length(name_of_test(test))))
  (α, Dn, Dcrit) = ks_parameters(test)
  p = p_value(test)
  result = p > test.α ? "FAIL TO REJECT H0" : "REJECT H0"
  pval = StatsBase.PValue(p)
  hypothesis(io, test)
  println(io, "Test report:")
  print(io, "  With $(Int((1 - test.α) * 100))% confidence: ", result)
  println(io)
  println(io, "  P-value:              ", pval)
  println(io, "Estimators:")
  print(io, "  D calculated:         ")
  show(io, test.Dn)
  println(io)
  print(io, "  D critical:           ")
  show(io, test.Dcrit)
  println(io)
  println(io, "Details:")
  show_params(io, test)
end