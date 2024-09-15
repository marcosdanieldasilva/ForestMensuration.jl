"""
Description:

ForestMensuration.jl is a Julia package that provides a comprehensive set of functions for performing dendrometric and forest inventory calculations. The package emphasizes ease of use, making it straightforward to:

- Conduct regressions: Linear models can be fitted with ease.
- Calculate tree and stand volume: Various methods are supported, including Huber, Smalian, and Newton.
- Perform forest inventories: Implementations of simple random sampling are available.
- Create frequency tables: Analyze the distribution of dendrometric variables.

Facilitates the analysis of dendrometric and forest data.
Performs complex calculations with simple commands.
Offers a user-friendly and intuitive interface.
"""
module ForestMensuration
  using
    CategoricalArrays,
    ColorTypes,
    DataFrames,
    Distributions,
    GLM,
    HypothesisTests,
    LinearAlgebra,
    PlotlyJS,
    Plots,
    Plots.PlotMeasures,
    RecipesPipeline,
    Reexport,
    StatsBase,
    StatsModels
    Tables

  import Plots: cgrad
  import StatsBase: fit, Histogram

  import StatsModels: asgn, missing_omit, Schema, TableRegressionModel

  @reexport using GLM

  include("structs-consts.jl")
  include("linear-regression.jl")
  include("multiple-linear-regression.jl")
  include("prediction.jl")
  include("regression-parameters.jl")
  include("criteria-functions.jl")
  include("plot_regression.jl")
  include("frequency-tables.jl")
  include("dendrometric-averages.jl")
  include("cubage.jl")
  include("inventory-report.jl")
  include("simple-casual-sampling.jl")
  include("site-classification.jl")
  include("show.jl")
  # include("graph-analysis.jl")

  export
    # Regression
    TableRegressionModel,
    regression,
    prediction,
    prediction!,
    criteria_table,
    criteria_selection,
    dendrometric_averages,
    plot_regression,
    #Cubage
    cubage,
    Smalian,
    Huber,
    Newton,
    # Frequency functions
    diametric_table,
    frequency_table,
    # Site classification
    hdom_classification,
    site_classification,
    site_table,
    # Forest Inventory
    simple_casual_sampling

end