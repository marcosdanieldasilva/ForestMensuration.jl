using ForestMensuration, CSV, DataFrames

data = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\AlturaDominante_vs_Idade.csv", DataFrame)

reg = regression(:Hd, :idade, data)


vol = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\volume_por_idade.csv", DataFrame)

reg2 = regression(:V, :DAP, :IDADE, vol, best=false);

ct =  criteria_table(reg2);

vol = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-volume-pinus.csv", DataFrame)


reg2 = regression(:vcc, :h, :d, vol, best=false)

ct =  criteria_table(reg2)



qreg = regression(:Hd, :idade, data, :parcelas)

s = site_classification(reg, 60)

site_table(reg, 60, 3)

site_table(reg, 60)

df = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-hipso.csv", DataFrame)

reg = regression(:ht, :dap, df)

p = graph(reg)

qreg = regression(:ht, :dap, df, :regiao)

qp = graph(qreg)

cub_data = CSV.read(
  "C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-cubagem-2.csv",
  DataFrame
)

cb = cubage(Smalian, :arvore, :h_d, :d_h, :casca, cub_data)


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