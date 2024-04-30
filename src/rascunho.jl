using ForestMensuration, CSV, DataFrames

df = CSV.read("C:\\Users\\marco\\OneDrive\\Documents\\Programação\\dados_CSV\\exemplo-hipso.csv", DataFrame)

reg = regression(:ht, :dap, df)

reg2 = regression(:ht, :dap, df, best=5)

cp = ForestMensuration._criteria_parameters(reg)


DataFrame(cp, ["RMSE", "Syx%", "Adj. R²", "Normality", "Homoscedasticity"])