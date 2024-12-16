import Base: show

# Custom show method for SiteAnalysis to display  the site_table and site_plot
function show(io::IO, analysis::SiteAnalysis)
  show(io, analysis.site_table)
  display(analysis.site_plot)
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

function show(io::IO, models::GroupedLinearModel)
  grouped_models = models.grouped_models
  isempty(grouped_models) && return show(io, grouped_models)
  summary(io, grouped_models), println(io, ":")
  for (k, v) in grouped_models
    println(io, "  $k => $v")
  end
end