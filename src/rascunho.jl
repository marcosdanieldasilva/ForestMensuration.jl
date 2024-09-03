using ForestMensuration

using CSV, DataFrames

df = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-hipso.csv", DataFrame)

reg = regression(:ht, :dap, df)

ct = criteria_table(reg)

cs = criteria_selection(reg)

plot_regression(cs)

reg2 = regression(:ht, :dap, df, :regiao)

ct2 = criteria_table(reg2)

cs2 = criteria_selection(reg2)


greg = regression(:ht, :dap, [:regiao], df)

greg2 = regression(:ht, :dap, [:regiao, :talhao], df)



greg2 = regression(:ht, :dap, :idade, [:regiao], df)

greg3 = regression(:ht, :dap, :idade, [:regiao], df)

greg4 = regression(:ht, :dap, :idade, [:regiao], df)

greg5 = regression(:ht, :dap, :idade, [:regiao], df, :talhao)


ct = criteria_table(reg)

reg2 = regression(:ht, :dap, df, model_type=LinearModel)

ct2 = criteria_table(reg2)

cs = criteria_selection(reg2)


function criteria_table(models::GroupedLinearModel, criteria::Symbol...; weight::Int=10)

  allowed_fields = [:adjr2, :syx, :rmse, :mae, :aic, :significance, :normality, :homoscedasticity]

  if isempty(criteria)
    selected_criteria = allowed_fields
  else
    selected_criteria = Any[s for s in criteria]
  end

  if !issubset(Set(selected_criteria), Set(allowed_fields))
    not_allowed = join(setdiff(selected_criteria, allowed_fields), ", :")
    allowed_msg = "\nAllowed fields are: :" * join(allowed_fields, ", :")
    throw(ArgumentError(":$not_allowed not allowed." * allowed_msg))
  end

  full_models = vcat(models, collect(models.grouped_models |> values))

  # # Generate the criteria parameters for each model in the input vector
  criteria_params = vcat(_criteria_parameters.(full_models)...)

  # Create a DataFrame from the full criteria parameters
  ct = DataFrame(criteria_params, allowed_fields)

  # # Filter the DataFrame columns based on the selected criteria
  ct = select(ct, selected_criteria)

  # # # Insert the model objects into the DataFrame
  insertcols!(ct, 1, "model" => full_models)
  insertcols!(ct, 1, join(models.group_names, " ") => vcat("Full Model", collect(keys(models.grouped_models))...))
  return ct
  # # Calculate ranks for each criterion
  # ranks = Dict()
  # if :adjr2 in selected_criteria
  #   ranks[:adjr2] = competerank(ct[!, "adjr2"], rev = true)
  # end
  # if :syx in selected_criteria
  #   ranks[:syx] = competerank(ct[!, "syx"])
  # end
  # if :rmse in selected_criteria
  #   ranks[:rmse] = competerank(ct[!, "rmse"])
  # end
  # if :mae in selected_criteria
  #   ranks[:mae] = competerank(ct[!, "mae"])
  # end
  # if :aic in selected_criteria
  #   ranks[:aic] = competerank(ct[!, "aic"])
  # end
  # if :significance in selected_criteria
  #   ranks[:significance] = competerank(ct[!, "significance"], rev = true) * weight
  # end
  # if :normality in selected_criteria
  #   ranks[:normality] = competerank(ct[!, "normality"], rev = true) * weight
  # end
  # if :homoscedasticity in selected_criteria
  #   ranks[:homoscedasticity] = competerank(ct[!, "homoscedasticity"], rev = true) * weight
  # end

  # # Combine ranks into a single score
  # combined_rank = sum([ranks[crit] for crit in selected_criteria])

  # # Insert the combined rank into the DataFrame
  # insertcols!(ct, 2, :rank => combined_rank)

  # # Sort the DataFrame by the combined rank
  # sort!(ct, :rank)

  # geral_models = [grouped_model.general_regression, grouped_model.qualy_regression]
  # geral_table = DataFrame(:name => ["General Model", "Qualy Model"], :model => geral_models)
  # geral_table = innerjoin(geral_table, criteria_table(geral_models, criteria..., best=false, weight=weight), on=:model)

  # grouped_models = vcat(models.grouped_models |> values |> collect)
  # group_table = DataFrame(:model => grouped_models)
  # insertcols!(group_table, 1, join(models.group_names, " ") => collect(keys(models.grouped_models)))
  # group_table = innerjoin(group_table, criteria_table(grouped_models, criteria..., best=false, weight=weight), on=:model)

  # return group_table
end                       

ct2 = criteria_table(reg)

cs2 = criteria_selection(reg)

reg2 = regression(:ht, :idade, :dap, df);

ct2 = criteria_table(reg2);

using StatsBase, Optim

using StatsModels, LinearAlgebra

function my_aic(model::StatsModels.TableRegressionModel)
  # Number of observations in the model
  n = nobs(model)
  # Degrees of freedom for residuals
  dof_resid = dof_residual(model)
  # The actual observed values (response variable)
  y = model.mf.data[1]
  # Predicted values from the model
  ŷ = prediction(model)
  # Residuals: the difference between observed and predicted values
  residual = y - ŷ
  # Deviance: sum of squared residuals
  devi = residual ⋅ residual
  # Log-likelihood of the model
  loglike = -n / 2 * (log(2π * devi / n) + 1)
  # Akaike Information Criterion (AIC): a measure of model quality
  AIC = -2 * loglike + 2 * dof_resid
  return AIC

end

best = ct[1, 1]

bic_glm(λ) = bic(glm(best.mf.f, df, Normal(), PowerLink(λ)));

optimal_bic = optimize(bic_glm, -1.0, 1.0);

round(optimal_bic.minimizer, digits = 5) # Optimal λ

best_opt = glm(best.mf.f, df, Normal(), PowerLink(optimal_bic.minimizer)) # Best model

round(optimal_bic.minimum, digits=5)


aic_glm(λ) = my_aic(glm(best.mf.f, df, Normal(), PowerLink(λ)));

optimal_aic = optimize(aic_glm, -1.0, 1.0);

round(optimal_aic.minimizer, digits = 5) # Optimal λ

best_opt2 = glm(best.mf.f, df, Normal(), PowerLink(optimal_aic.minimizer)) # Best model

round(optimal_aic.minimum, digits=5)

ct3 = criteria_table([best, best_opt, best_opt2])

data = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\AlturaDominante_vs_Idade.csv", DataFrame)

reg = regression(:Hd, :idade, data)

ct = criteria_table(reg)


using StatsModels

d, ht = term.([:dap, :ht])

tems = ForestMensuration.generate_combined_terms(d, ht)

data = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\AlturaDominante_vs_Idade.csv", DataFrame)

reg = regression(:Hd, :idade, data, model_type=LinearModel)

ct = criteria_table(reg)

best = site_table(reg[1], 60, 2)



reg2 = regression(:Hd, :idade, data)

ct2 = criteria_table(reg2)

best2 = site_table(reg2[1], 60, 2)


reg = regression(:Hd, :idade, data)

ct =  criteria_table(reg);

reg2 = regression(:Hd, :idade, data, model_type=LinearModel)

ct2 =  criteria_table(reg2);

best = criteria_table(reg)[end, 1]

vol = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\volume_por_idade.csv", DataFrame)

reg = regression(:V, :DAP, vol)

ct =  criteria_table(reg);


reg2 = regression(:V, :DAP, vol, model_type=LinearModel)

ct2 =  criteria_table(reg2);

reg2 = regression(:V, :DAP, :IDADE, vol)

ct2 =  criteria_table(reg2);

reg3 = regression(:V, :DAP, :IDADE, vol)

ct3 =  criteria_table(reg3);


reg4 = regression(:V, :DAP, :IDADE, vol)

ct4 =  criteria_table(reg4);


best = reg2[2]

X = best.mm.m;
Y =  best.model.rr.y;


fitted_model = GLM.fit(GeneralizedLinearModel, X, Y, Normal(), LogLink())

fitted_model2 = GLM.fit(GeneralizedLinearModel, X, Y, Normal(), IdentityLink())

using GLM

reg_glm = glm(best.mf.f, vol, Normal(), LogLink())

reg3 = regression(:V, :DAP, :IDADE, vol, :ARV, best=false);

ct3 =  criteria_table(reg3);

vol = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-volume-pinus.csv", DataFrame)

reg = regression(:vcc, :h, :d, vol, best=5);

ct =  criteria_table(reg); 

reg2 = regression(:vcc, :h, :d, vol, best=false);

ct2 =  criteria_table(reg2);


# julia> ct
# 5×7 DataFrame
#  Row │ Adj. Model                         Rank   RMSE     Syx%     Adj. R²  Normality  Homoscedasticity 
#      │ FittedLin…                         Int64  Float64  Float64  Float64  Float64    Float64
# ─────┼──────────────────────────────────────────────────────────────────────────────────────────────────
#    1 │ log(vcc) = -3.5475 + 0.01839 * l…     23   0.0628    14.91   0.7178        1.0               0.0
#    2 │ log(vcc) = -3.5475 - 0.0001032 *…     23   0.0628    14.91   0.7178        1.0               0.0
#    3 │ log(vcc) = -5.1006 - 0.002613 * …     29   0.0631    14.98   0.7154        1.0               0.0
#    4 │ log(vcc) = -5.1006 + 0.1204 * h …     29   0.0631    14.98   0.7154        1.0               0.0
#  5 │ vcc = -0.2955 + 0.001035 * log(h…     35   0.0632    15.14   0.7092        1.0               0.0

# Row │ Adj. Model                         Rank   RMSE     Syx%     Adj. R²  Normality  Homoscedasticity 
#      │ FittedLin…                         Int64  Float64  Float64  Float64  Float64    Float64
# ─────┼──────────────────────────────────────────────────────────────────────────────────────────────────
#    1 │ vcc = -9237.4158 + 6478.0 * log(…     61   0.0554    13.9    0.7548        1.0               1.0
#    2 │ vcc = -3654.455 + 2568.0 * log(h…     95   0.0555    13.93   0.7538        1.0               1.0
#    3 │ vcc = -1366.8106 + 962.3 * log(h…     99   0.0555    13.93   0.7537        1.0               1.0
#    4 │ vcc = -318.092 - 10.19 * h + 173…    105   0.0555    13.94   0.7535        1.0               1.0
#    5 │ vcc = -421.1801 - 13.42 * h + 22…    105   0.0555    13.94   0.7535        1.0               1.0
#    6 │ vcc = -902.1231 - 28.61 * h + 49…    107   0.0555    13.94   0.7534        1.0               1.0

reg3 = regression(:vcc, :h, :d, vol, best=false);

ct3 =  criteria_table(reg3);

reg2 = regression(:vcc, :h, :d, vol, best=false);

ct =  criteria_table(reg2);

qreg = regression(:Hd, :idade, data, :parcelas)

s = site_classification(reg, 60)

site_table(reg, 60, 3)

site_table(reg, 60)



ct2 = criteria_table(reg2)

reg3 = regression(:ht, :dap, df);

ct3 = criteria_table(reg3)

df2 = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-inventario.csv", DataFrame);


reg2 = regression(:ht, :dap, df2);

ct2 = criteria_table(reg2);

cs2 = criteria_selection(reg2)

reg2 = regression(:ht, :dap, df2, model_type=Tables)

ct2 = criteria_table(reg2);

cs2 = criteria_selection(reg2)

prediction!(cs2, df2)

sort!(df2, :dap);
scatter(df2.dap, df2.ht);
plot!(df2.dap, df2.ht_predict)

best = ct2[1, 1]

sitio = regression(:hdom, :idade, :dap, d3);

p = graph(reg)

qreg = regression(:ht, :dap, df, :regiao, :parcela)

qp = graph(qreg)

cub_data = CSV.read(
  "C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-cubagem-2.csv",
  DataFrame
)

cb = cubage(Smalian, :arvore, :h_d, :d_h, :casca, cub_data, 6)


gdf = groupby(cub_data, :arvore)

ForestMensuration._diameter_interpolation(gdf[1].h_d[end] * 0.1, gdf[1].h_d*1, gdf[1].d_h*1)


# julia> @btime regression(:ht, :dap, df)
# 9.203 ms (14741 allocations: 19.69 MiB)
# 1-element Vector{FittedLinearModel}:
#  log(ht - 1.3) = 1.2946 + 0.1516 * dap - 0.002818 * dap ^ 2
reg = regression(:ht, :dap, df)


reg2 = regression(:ht, :dap, df, best=5)

cp = ForestMensuration._criteria_parameters(reg)


DataFrame(cp, ["RMSE", "Syx%", "Adj. R²", "Normality", "Homoscedasticity"])

V = [1	11.0	220
2	10.0	200
3	10.5	210
4	12.0	240
5	9.8	196
6	12.3	246
7	11.9	238
8	12.0	240
9	13.2	264
10	9.6	192
11	9.9	198
12	10.3	206]

v = DataFrame(V, [:parcela, :m³, :volume_ha])

ass = simple_casual_sampling(v.volume_ha, 1, 550)


@btime regression(:ht, :dap, df, :regiao)

regression(:ht, :dap, df, :regiao)

qreg = @btime regression(:ht, :dap, df, :regiao, :talhao);

qreg2 = @btime regression(:ht, :dap, df, :regiao, effect=:interactive);