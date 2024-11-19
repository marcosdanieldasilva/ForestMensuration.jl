"""
    const MixTerm = Union{AbstractTerm,Tuple{AbstractTerm,Vararg{AbstractTerm}}}

Union type representing a mixed term, which can be a single `AbstractTerm` or a tuple of `AbstractTerm`s.
"""
const MixTerm = Union{AbstractTerm,Tuple{AbstractTerm,Vararg{AbstractTerm}}}

"""
    const β0 = InterceptTerm{true}()

Represents an intercept term for linear models.
"""
const β0 = InterceptTerm{true}()

"""
    const emptySchema = Schema()

Represents a default schema used for modeling operations.
"""
const emptySchema = Schema()

"""
    abstract type CubingMethod
  
Abstract type representing a method for cubing (calculating volume).
  # Subtypes
  - Smalian
  - Huber
  - Newton

"""
abstract type CubingMethod end

"""
    abstract type Smalian <: CubingMethod

Smalian Method:
  The Smalian method measures diameters or circumferences at the ends of each section and calculates 
  the total volume by:
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
    abstract type Huber <: CubingMethod

Huber Method:
  The Huber method measures the diameter or circumference at the midpoint of the section, and the volume
   is determined by:
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
    abstract type Newton <: CubingMethod

Newton Method:
  The Newton method involves measuring at 3 positions along each section (at the ends and in the middle 
  of the logs). Therefore, it is a more laborious method than the others, but the estimated volume will 
  be more accurate.
  
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
    struct ModelEquation

Define ModelEquation struct to store the regression results
  # Fields
  - output::String
  - model::TableRegressionModel
"""
struct ModelEquation
  output::String
  model::TableRegressionModel
  # Inner constructor to initialize `output` based on `model`
  function ModelEquation(model::TableRegressionModel)
    # Get coefficients and terms from the model
    β = coef(model)
    n = length(β)
    output = string(StatsModels.coefnames(model.mf.f.lhs)) * " = $(round(β[1], digits = 6))"

    for i in 2:n
      term = coefnames(model)[i]
      product = string(round(abs(β[i]), sigdigits=6)) * " * " * term
      output *= signbit(β[i]) ? " - $(product)" : " + $(product)"
    end

    # Return the constructed RegressionEquation object
    new(output, model)
  end
end

"""
    struct SiteAnalysis

Define SiteAnalysis struct to store the analysis results
  # Fields
  - site_table::DataFrame
  - site_plot::Plots.Plot
"""
struct SiteAnalysis
  site_table::DataFrame
  site_plot::Plots.Plot
end