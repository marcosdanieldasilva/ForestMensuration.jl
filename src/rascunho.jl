using ForestMensuration, CSV, DataFrames

df = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-hipso.csv", DataFrame)

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

