"""
Description:

ForestMensuration.jl is a Julia package that provides a comprehensive set of functions for performing dendrometric/biometric and forest inventory calculations. The package emphasizes ease of use, making it straightforward to:

Conduct regressions: Linear, nonlinear, and mixed-effects models can be fitted with ease.
Calculate tree and stand volume: Various methods are supported, including Huber, Smalian, and Newton.
Perform forest inventories: Implementations of simple random and stratified sampling are available.
Create frequency tables: Analyze the distribution of dendrometric variables.
Work with probability distributions: Calculate PDF and CDF for various distributions.

Facilitates the analysis of dendrometric and forest data.
Performs complex calculations with simple commands.
Offers a user-friendly and intuitive interface.
"""
module ForestMensuration
  using DataFrames, Distributions, LinearAlgebra, HypothesisTests, StatsBase, StatsModels, Tables
  import StatsBase: dof_residual, dof, nobs, aicc, aic, bic, coef

  
  include("structs-consts.jl")
  include("regression-variables.jl")
  include("linear-regression.jl")
  include("regression-parameters.jl")
  include("frequency-tables.jl")
  include("cubage.jl")
  include("inventory-report.jl")
  include("simple-casual-sampling.jl")
  include("show.jl")

  export
    # regression structures
    FittedLinearModel,
    # cubage methods
    Smalian,
    Huber,
    Newton,
    # funtions
    adjr2,
    aic,
    aicc,
    bark_factor,
    bic,
    coef,
    coef_table,
    confidence_interval,
    criteria_table,
    cubage,
    deviance,
    diametric_table,
    dispersion,
    dof,
    dof_residual,
    fit_regression,
    frequency_table,
    homoscedasticity,
    loglikelihood,
    modelmatrix,
    n_coef,
    nobs,
    normality,
    nulldeviance,
    predict,
    r2,
    regression,
    residuals,
    stderror,
    syx,
    syx_in_percentage,
    simple_casual_sampling,
    stratified_sampling

end