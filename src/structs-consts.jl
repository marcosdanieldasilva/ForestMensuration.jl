"""
Represents a fitted linear model.

# Fields
- `formula::F`: The formula used for the model.
- `data::D`: The data frame containing the data.
- `β::Array{T, 1}`: The regression coefficients.
- `σ²::T`: The variance of residuals.
- `RMSE::T`: The root mean squared error.
- `chol::C`: The Cholesky decomposition of X'X.
"""
struct FittedLinearModel{F <: FormulaTerm, D <: AbstractDataFrame, T <: Float64, C <: Cholesky{Float64, Matrix{Float64}}}
  formula::F
  data::D
  β::Array{T, 1}
  σ²::T
  RMSE::T
  chol::C
end

"""
Creates a new `FittedLinearModel` instance.

# Arguments
- `formula::F`: The formula used for the model.
- `data::D`: The data frame containing the data.
- `β::Array{T, 1}`: The regression coefficients.
- `σ²::T`: The variance of residuals.
- `RMSE::T`: The root mean squared error.
- `chol::C`: The Cholesky decomposition of X'X.

# Returns
- `FittedLinearModel{F, D, T, C}`: A new fitted linear model.
"""
function FittedLinearModel(formula::F, data::D, β::Array{T, 1}, σ²::T, RMSE::T, chol::C) where {F, D, T, C}
  return FittedLinearModel{F, D, T, C}(formula, data, β, σ², RMSE, chol)
end

"""
Union type representing a mixed term, which can be a single `AbstractTerm` or a tuple of `AbstractTerm`s.
"""
const MixTerm = Union{AbstractTerm, Tuple{AbstractTerm, Vararg{AbstractTerm}}}

"""
Represents an intercept term for linear models.
"""
const β0 = InterceptTerm{true}()

"""
Abstract type representing a method for cubing (calculating volume).
"""
abstract type CubingMethod end

"""
Smalian Method:
  The Smalian method measures diameters or circumferences at the ends of each section and calculates the total volume by:
  - Vt = v0 + Σi=1:n(vi) + vt
  - v0 = g0 * l0
  - vi = (gi + gi+1)/2 * li
  - vt = (1/3) * gn * ln
  Where:
  - v0 = volume of the stump;
  - vi = volume of intermediate sections;
  - vt = volume of the cone;
  - g = basal area;
  - l = length.
"""
abstract type Smalian <: CubingMethod end

"""
Huber Method:
  The Huber method measures the diameter or circumference at the midpoint of the section, and the volume is determined by:
  - v = v0 + Σi=1:n(vi) + vt
  - vi = gi * li
  Where:
  - v0 = volume of the stump;
  - vi = volume of intermediate sections;
  - vt = volume of the cone;
  - g = basal area;
  - l = length.
"""
abstract type Huber <: CubingMethod end

"""
Newton Method:
  The Newton method involves measuring at 3 positions along each section (at the ends and in the middle of the logs). Therefore, it is a more laborious method than the others, but the estimated volume will be more accurate.
  
  - v = v0 + Σi=1:n(vi) + vt
  - vi = (gi + gm + gi+1)/2 * li
  Where:
  - v0 = volume of the stump;
  - vi = volume of intermediate sections;
  - vt = volume of the cone;
  - g = basal area;
  - gm = basal area at the midpoint of the section;
  - l = length.
"""
abstract type Newton <: CubingMethod end

"""
Represents a report for stratified sampling in forest inventory.

# Fields
- `anova::AbstractDataFrame`: ANOVA table.
- `auxiliary_table::AbstractDataFrame`: Auxiliary table with additional statistics.
- `result_table::AbstractDataFrame`: Result table with final estimates.
"""
struct StratifiedReport
  anova::AbstractDataFrame
  auxiliary_table::AbstractDataFrame
  result_table::AbstractDataFrame
end

# Define SiteAnalysis struct to store the analysis results
struct SiteAnalysis
  site_table::DataFrame
  site_plot::Plots.Plot
end