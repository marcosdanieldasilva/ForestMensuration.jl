"""
Generate an inventory report.

# Arguments
- `x̅::F`: Mean volume per plot.
- `cv::F`: Coefficient of variation.
- `Ttab::F`: T-value for the confidence level.
- `s²x̅::F`: Mean variance.
- `sx̅::F`: Standard error.
- `plot_area::Union{I, F, S}`: Area of the plot.
- `f::F`: Sampling fraction.
- `population::String`: Population type ("infinite" or "finite").
- `n::Union{I, T}`: Number of measured plots.
- `required_plots::Union{I, T}`: Number of required plots.
- `missing_plots::Union{I, T}`: Number of missing plots.
- `N::I`: Total number of plots.
- `lg::Symbol`: Language for the report (":pt" for Portuguese, ":en" for English).

# Returns
- `DataFrame`: The inventory report.
"""
function _inventory_report(x̅::F, cv::F, Ttab::F, s²x̅::F, sx̅::F, plot_area::Union{I, F, S}, f::F, population::String, n::Union{I, T}, required_plots::Union{I, T}, missing_plots::Union{I, T}, N::I, lg::Symbol) :: DataFrame where F <: Float64 where I <: Int where S <: String where T <: Tuple
  hectare_volume = plot_area == "u.a." ? NaN64 : round(x̅ / plot_area, digits = 3)
  absolute_error = Ttab * sx̅
  relative_error = round(absolute_error / x̅ * 100, digits = 2)
  total_volume = round(x̅ * N, digits = 3)
  ic_lower = round(total_volume - N * absolute_error, digits = 3)
  ic_upper = round(total_volume + N * absolute_error, digits = 3)
  hectare_volume === NaN64 ? N = Inf : nothing
  unit = "m³/" * string(plot_area) * "ha"
  
  values = [
    x̅, round(cv, digits = 2), s²x̅, sx̅, absolute_error, 
    relative_error, hectare_volume, total_volume, ic_lower, ic_upper, 
    f, n, required_plots, missing_plots, N
  ]
  units = [
    unit, "%", "($(unit))²", unit, unit, 
    "%", "m³ha⁻¹", "m³", "m³", "m³", 
    population, "n", "n", "n", "N"
  ]
  
  if lg == :pt
    return DataFrame(
      Parâmetros = [
        "volume por parcela", "coeficiente de variação", "variância média", "erro padrão", "erro absoluto",
        "erro relativo", "volume por hectare", "Volume Total", "intervalo de confiança inferior", "intervalo de confiança superior",
        "população", "parcelas medidas", "parcelas necessárias", "parcelas faltantes", "parcelas possíveis"
      ],
      Valores = values,
      Unidades = units
    )
  else
    return DataFrame(
      Parameters = [
        "plot volume", "coefficient of variation", "mean variance", "standard error", "absolute error",
        "relative error", "hectare volume", "Total Volume", "confidence interval lower", "confidence interval upper",
        "population", "measured plots", "required plots", "missing plots", "possible plots"
      ],
      Values = values,
      Units = units
    )
  end
end