using ForestMensuration, CSV, DataFrames

using StatsModels

d, ht = term.([:dap, :ht])

tems = ForestMensuration.generate_combined_terms(d, ht)

data = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\AlturaDominante_vs_Idade.csv", DataFrame)

reg = regression(:Hd, :idade, data)


vol = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\volume_por_idade.csv", DataFrame)

reg2 = regression(:V, :DAP, :IDADE, vol, best=false);

ct2 =  criteria_table(reg2);

reg2 = regression(:IDADE, :DAP, :V, vol);


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

df = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-hipso.csv", DataFrame)

reg = regression(:ht, :dap, df,best=false);

sitio = regression(:hdom, :idade, :dap, df);

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