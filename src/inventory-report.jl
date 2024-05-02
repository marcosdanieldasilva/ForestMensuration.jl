function _inventory_report(x̅::F, cv::F, Ttab::F, s²x̅::F, sx̅::F, plot_area::Union{I, F, S}, f::F, population::String, n::Union{I, T}, required_plots::Union{I, T}, missing_plots::Union{I, T}, N::I, lg::Symbol) ::DataFrame where F <: Float64 where I <: Int where S <: String where T <: Tuple
  plot_area == "u.a." ? hectare_volume = NaN64 : hectare_volume = round(x̅ / plot_area, digits = 3)
  absolute_error = Ttab * sx̅
  relative_error = round(absolute_error / x̅ * 100, digits = 2)
  total_volume = round(x̅ * N, digits = 3)
  ic_lower = round(total_volume - N * absolute_error, digits = 3)
  ic_upper = round(total_volume + N * absolute_error, digits = 3)
  hectare_volume === NaN64 ? N = Inf : nothing
  unit = "m³/"*"$plot_area"*"ha"
  V = [
    x̅, round(cv, digits = 2), s²x̅, sx̅, absolute_error, 
    relative_error, hectare_volume, total_volume, ic_lower, ic_upper, 
    f, n, required_plots, missing_plots, N
  ]
  U = [
    unit, "%", "($(unit))²", unit, unit, 
    "%", "m³ha⁻¹", "m³", "m³", "m³", 
    population, "n", "n", "n", "N"
  ]
  return if lg == :pt
    DataFrame(
      Parâmetros = [
        "volume por parcela", "coeficiente de variação", "variância média", "erro padrão", "erro absoluto",
        "erro relativo", "volume por hectare", "Volume Total", "intervalo de confiança inferior", "intervalo de confiança superior",
        "população", "parcelas medidas", "parcelas necessárias", "parcelas faltantes", "parcelas possíveis"
      ],
      Valores = V,
      Unidades = U
    )
  else
    DataFrame(
      Parameters = [
      "plot volume", "coefficient of variation", "mean variance", "standard error", "absolute error",
      "relative error", "hectare volume", "Total Volume", "confidence interval lower", "confidence interval upper",
      "population", "measured plots", "required plots", "missing plots", "possible plots"
      ],
      Values = V,
      Units = U
    )
  end
end