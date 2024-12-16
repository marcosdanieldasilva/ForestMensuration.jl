"""
Description:

ForestMensuration.jl is a Julia package that provides a comprehensive set of functions for performing dendrometric and forest inventory calculations. Designed with ease of use in mind, the package offers tools that simplify complex forestry calculations, making it straightforward to:

- **Conduct regressions**: Easily fit linear models, including transformations and qualitative variables, to find the best relationships in your forestry data.
- **Calculate tree and stand volume (cubage)**: Support for various methods such as Huber, Smalian, and Newton allows precise calculation of tree and stand volumes.
- **Classify site productivity**: Implement site classification methods to assess the quality and productivity of forest sites.
- **Compute dendrometric averages**: Calculate essential dendrometric metrics like mean diameter, quadratic mean diameter, and others to understand stand structure.
- **Create frequency tables**: Generate frequency and diametric tables to analyze the distribution of dendrometric variables such as diameter and height.
- **Perform forest inventory**: Implement simple random sampling methods to estimate forest parameters and support forest management decisions.

The package facilitates the analysis of dendrometric and forest data, performs complex calculations with simple commands, and offers a user-friendly and intuitive interface.
"""
module ForestMensuration
using
  CategoricalArrays,
  ColorTypes,
  DataFrames,
  Distributions,
  HypothesisTests,
  LinearAlgebra,
  PlotlyJS,
  Plots,
  Plots.PlotMeasures,
  RecipesPipeline,
  Reexport,
  ScientificTypes,
  StatsBase,
  StatsModels,
  Tables

import Plots: cgrad
import StatsBase: fit, Histogram

include("structs_consts.jl")
include("goodness_of_fit_test.jl")
include("linear_regression.jl")
include("prediction.jl")
include("regression_parameters.jl")
include("criteria_functions.jl")
include("plot_regression.jl")
include("frequency_tables.jl")
include("dendrometric_averages.jl")
include("cubage.jl")
include("inventory_report.jl")
include("simple_casual_sampling.jl")
include("site_classification.jl")
include("show.jl")

export
  # Regression
  LinearModel,
  GroupedLinearModel,
  regression,
  predict,
  predict!,
  coef_table,
  criteria_table,
  criteria_selection,
  plot_regression,
  #Cubage
  artificial_form_factor,
  bark_factor,
  bole_volume,
  cone_volume,
  cubage,
  cylinder_volume,
  natural_form_factor,
  quotient_form,
  Smalian,
  Huber,
  Newton,
  # Frequency and Statistic functions
  basal_area,
  dendrometric_averages,
  diametric_table,
  frequency_table,
  # Site classification
  hdom_classification,
  site_classification,
  site_table,
  # Forest Inventory
  simple_casual_sampling

end