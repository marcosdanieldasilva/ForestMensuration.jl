function _infinite_proportional_allocation(ps²::Real, ni::Real, e::Real, α::Real)
  Ttab1 = -quantile(TDist(ni - 1), (1 - α) / 2)
  iter1 = (Ttab1^2 * ps²) / e^2
  Ttab2 = -quantile(TDist(iter1 - 1), (1 - α) / 2)
  iter2 = (Ttab2^2 * ps²) / e^2
  if round(ni, digits=3) != round(iter1, digits=3) && round(ni, digits=3) != round(iter2, digits=3)
    _infinite_proportional_allocation(ps², iter1, e, α)
  else
    return ceil(Integer, ni)
  end
end

function _finite_proportional_allocation(ps²::Real, ni::Real, N::Real, e::Real, α::Real)
  Ttab1 = -quantile(TDist(ni - 1), (1 - α) / 2)
  iter1 = (Ttab1^2 * ps²) / (e^2 + Ttab1^2 * (ps² / N))
  Ttab2 = -quantile(TDist(iter1 - 1), (1 - α) / 2)
  iter2 = (Ttab2^2 * ps²) / (e^2 + Ttab2^2 * (ps² / N))
  if round(ni, digits=3) != round(iter1, digits=3) && round(ni, digits=3) != round(iter2, digits=3)
    _finite_proportional_allocation(ps², iter1, N, e, α)
  else
    return ceil(Integer, ni)
  end
end

_diffn(t::NTuple{N,T}) where {N,T} = ntuple(i -> t[i] - t[i+1], N - 1)

_diff(t::NTuple{N,T}) where {N,T} = ntuple(i -> t[i+1] - t[i], N - 1)

function _anova(h::Symbol, q::Symbol, data::AbstractDataFrame)
  new_data = dropmissing(data[:, [h, q]])
  formula = (term(h) ~ term(1) + term(q)) |> (formula -> apply_schema(formula, StatsModels.schema(formula, new_data, Dict(h => ContinuousTerm, q => CategoricalTerm))))
  Y, X = modelcols(formula, new_data)
  SST = sum(abs2.(new_data[!, 1] .- mean(new_data[!, 1])))
  n = size(new_data, 1)
  ŷ = similar(Vector{Float64}, n)
  β = X'Y
  # Compute the Cholesky decomposition of X'X for optimization
  chol = cholesky!(X'X)
  # Calculate the coefficients of the fitted regression
  ldiv!(chol, β)
  # Calculate the predicted values
  mul!(ŷ, X, β)
  # Calculate the residuals values
  residual = Y - ŷ
  # Sum Of Square Error
  SSR = residual ⋅ residual
  dof = (2, length(β) + 1)
  Δdof = _diff(dof)
  dof_residual = @. Int(n - dof + 1)
  MSR1 = _diffn((SST, SSR)) ./ Δdof
  MSR2 = ((SST, SSR) ./ dof_residual)
  forward = dof[1] <= dof[2]

  if forward
    MSR2 = MSR2[2:end]
    dfr_big = dof_residual[2:end]
  else
    MSR2 = MSR2[1:end-1]
    dfr_big = dfr[1:end-1]
  end

  fstat = (NaN, (MSR1 ./ MSR2)...)
  pval = (NaN, ccdf.(FDist.(abs.(Δdof), dfr_big), abs.(fstat[2:end]))...)
  DataFrame("DOF" => [dof...], "ΔDOF" => ["", Δdof[1]], "SSR" => [SST, SSR], "F*" => [fstat...], "p(>F)" => [pval...])
end

function _auxiliary_table(stratum::Symbol, volume::Symbol, N::Vector, data::AbstractDataFrame)
  @chain data begin
    groupby(_, stratum)
    combine(_) do df
      (n=length(df[!, stratum]), x̅=mean(df[!, volume]), s²=var(df[!, volume]), s=std(df[!, volume]))
    end
    insertcols!(_, 6, :p => N ./ sum(N))
    insertcols!(_, 7, :ps => _.p .* _.s)
    insertcols!(_, 8, :ps² => _.p .* _.s²)
  end
end

function _difference_and_remainders(measured::Tuple, ideal::Tuple)::Tuple
  difference = ideal .- measured
  result = max.(difference, 0)
  return result
end

function stratified_sampling(stratum::Symbol, volume::Symbol, plot_area::Real, total_area::Vector, data; e::Real=10, α::Real=0.95, lg::Symbol=:pt)

  anova = _anova(volume, stratum, data)
  N = round(Integer, sum(total_area) / plot_area)
  auxiliary_table = _auxiliary_table(stratum, volume, total_area, data)
  x̅st = sum(auxiliary_table.p .* auxiliary_table.x̅)
  s²x̅st = sum(auxiliary_table.p .^ 2 .* (auxiliary_table.s² / auxiliary_table.n) .* (1 .- (auxiliary_table.n ./ total_area)))
  sx̅st = √s²x̅st
  cv = sx̅st / x̅st * 100
  E = (e / 100) * x̅st
  n = tuple(auxiliary_table.n...)

  required_plots = tuple(round.(Integer, _finite_proportional_allocation(sum(auxiliary_table.ps²), sum(auxiliary_table.n), N, E, α) .* auxiliary_table.p)...)
  missing_plots = _difference_and_remainders(n, required_plots)

  Ttab = begin
    gh = (total_area .* (total_area .- auxiliary_table.n)) ./ auxiliary_table.n
    en = ceil(Integer, sum(gh .* auxiliary_table.s²)^2 / sum((gh .^ 2 .* auxiliary_table.s² .^ 2) ./ (auxiliary_table.n .- 1)))
    -quantile(TDist(en), (1 - α) / 2)
  end

  f = 1 - (sum(auxiliary_table.n) / (sum(total_area) / plot_area))
  population = lg == :pt ? "infinita" : "infinite"

  result_table = ForestMensuration._inventory_report(x̅st, cv, Ttab, s²x̅st, sx̅st, plot_area, f, population, n, required_plots, missing_plots, N, lg)

  StratifiedReport(anova, auxiliary_table, result_table)
end