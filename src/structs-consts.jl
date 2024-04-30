struct FittedLinearModel{F <: FormulaTerm, D <: AbstractDataFrame, T <: Float64, C <: Cholesky{Float64, Matrix{Float64}}}
  formula::F
  data::D
  β::Array{T, 1}
  σ²::T
  RMSE::T
  chol::C
end

function FittedLinearModel(formula::F, data::D, β::Array{T, 1}, σ²::T, RMSE::T, chol::C) where {F, D, T, C}
  return FittedLinearModel{F, D, T, C}(formula, data, β, σ², RMSE, chol)
end

const MixTerm = Union{AbstractTerm, Tuple{AbstractTerm, Vararg{AbstractTerm}}}

const β0 = InterceptTerm{true}()

abstract type CubingMethod end

"""
Smalian Method:
  The Smalian method measures diameters or circumferences at the ends of each section and calculates the total volume by:
  - Vt = v0 + Σi=1:n(vi) + vt
  - v0 = g0 * l0
  - vi = (gi+gi+1)/2 * li
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

struct StratifiedReport
  anova::AbstractDataFrame
  auxiliary_table::AbstractDataFrame
  result_table::AbstractDataFrame
end