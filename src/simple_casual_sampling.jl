# Calculates the sample size for an infinite population.
# Arguments
# - `t::Real`: T-value for the given confidence level.
# - `cv::Real`: Coefficient of variation.
# - `e::Real`: Desired error margin.
# Returns
# - `Float64`: The sample size for an infinite population.
_infinite_sample(t::Real, cv::Real, e::Real) = (t^2 * cv^2) / e^2

# Calculates the sample size for a finite population.
# Arguments
# - `t::Real`: T-value for the given confidence level.
# - `cv::Real`: Coefficient of variation.
# - `e::Real`: Desired error margin.
# - `N::Int64`: Population size.
# Returns
# - `Float64`: The sample size for a finite population.
_finite_sample(t::Real, cv::Real, e::Real, N::Int64) = (t^2 * cv^2) / (e^2 + ((t^2 * cv^2) / N))

# Calculates the required sample size for an infinite population using an iterative approach with a 
#maximum iteration limit to prevent infinite loops.
# Arguments
# - `cv::Real`: Coefficient of variation.
# - `ni::Real`: Initial sample size.
# - `e::Real`: Desired error margin.
# - `α::Real`: Confidence level.
# - `max_iterations::Int=1000`: Maximum number of iterations (default is 1000).
# Returns
# - `Int64`: The required sample size for an infinite population.
function _infinite_population(cv::Real, ni::Real, e::Real, α::Real; max_iterations::Int=1000)
  iteration = 0
  while iteration < max_iterations
    Ttab1 = -quantile(TDist(ni - 1), (1 - α) / 2)
    iter1 = _infinite_sample(Ttab1, cv, e)
    Ttab2 = -quantile(TDist(iter1 - 1), (1 - α) / 2)
    iter2 = _infinite_sample(Ttab2, cv, e)
    if round(ni, digits=3) == round(iter1, digits=3) || round(ni, digits=3) == round(iter2, digits=3)
      return ceil(Int, ni)
    end
    ni = iter1
    iteration += 1
  end
  error("Maximum number of iterations ($max_iterations) reached without convergence.")
end

# Calculates the required sample size for a finite population using an iterative approach with a maximum
# iteration limit to prevent infinite loops.
# Arguments
# - `cv::Real`: Coefficient of variation.
# - `ni::Real`: Initial sample size.
# - `N::Int`: Population size.
# - `e::Real`: Desired error margin.
# - `α::Real`: Confidence level.
# - `max_iterations::Int=1000`: Maximum number of iterations (default is 1000).
# Returns
# - `Int64`: The required sample size for a finite population.
function _finite_population(cv::Real, ni::Real, N::Int, e::Real, α::Real; max_iterations::Int=1000)
  iteration = 0
  while iteration < max_iterations
    Ttab1 = -quantile(TDist(ni - 1), (1 - α) / 2)
    iter1 = _finite_sample(Ttab1, cv, e, N)
    Ttab2 = -quantile(TDist(iter1 - 1), (1 - α) / 2)
    iter2 = _finite_sample(Ttab2, cv, e, N)
    if round(ni, digits=3) == round(iter1, digits=3) || round(ni, digits=3) == round(iter2, digits=3)
      return ceil(Int, ni)
    end
    ni = iter1
    iteration += 1
  end
  error("Maximum number of iterations ($max_iterations) reached without convergence.")
end


"""
    simple_casual_sampling(volume::Vector{<:Real}, plot_area::Real, total_area::Real; 
                           e::Real=10, α::Real=0.95, lg::Symbol=:en)

Performs simple random sampling for forest inventory analysis with specified plot area and total area.

# Description:

The `simple_casual_sampling` function calculates various statistical parameters for a forest inventory using simple random sampling methodology. It is designed to estimate the total volume of timber within a forest area by analyzing volume measurements from sample plots.

Simple random sampling is a statistical method where each unit (in this case, a plot) has an equal chance of being selected. This method is commonly used in forest inventories to estimate parameters like mean volume per hectare, total volume, and confidence intervals.

The function computes key statistics such as:

- Mean plot volume
- Coefficient of variation
- Variance and standard error of the mean
- Absolute and relative errors
- Volume per hectare
- Total estimated volume
- Confidence intervals for the total volume
- Required number of plots to achieve a desired error margin
- Number of missing plots to be measured

# Arguments:

- `volume::Vector{<:Real}`: A vector containing the volume measurements from each sampled plot (e.g., in cubic meters per plot).
- `plot_area::Real`: The area of each sample plot (e.g., in hectares). All plots are assumed to have the same area.
- `total_area::Real`: The total area of the forest or stand being inventoried (e.g., in hectares).
- `e::Real=10`: The desired relative error margin as a percentage (default is 10%). This represents the maximum acceptable error in the estimate relative to the true mean.
- `α::Real=0.95`: The confidence level for the statistical estimates (default is 95%). This determines the confidence interval for the total volume estimate.
- `lg::Symbol=:en`: Language for the report output (default is English). Supported options are `:en` for English and `:pt` for Portuguese.

# Returns:

- `report::DataFrame`: A DataFrame containing the inventory report with calculated statistical parameters and their values, including units of measurement.

# Example:

```Julia
# Define the volumes measured in each sample plot (e.g., in cubic meters)
julia> v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1];

# Define the area of each plot (in hectares)
julia> plot_area = 0.05 # hectares

# Define the total area of the forest (in hectares)
julia> total_area  = 10 # hectares

# Perform the simple random sampling analysis
julia> simple_casual_sampling(v, plot_area, total_area; e=10, α=0.95, lg=:en)
15×3 DataFrame
 Row │ Parameters                 Values      Units        
     │ String                     Float64     String       
─────┼─────────────────────────────────────────────────────
   1 │ plot volume                  441.691   m³/0.05ha
   2 │ coefficient of variation      10.03    %
   3 │ mean variance                168.498   (m³/0.05ha)²
   4 │ standard error                12.9807  m³/0.05ha
   5 │ absolute error                28.9227  m³/0.05ha
   6 │ relative error                 6.55    %
   7 │ hectare volume              8833.82    m³ha⁻¹
   8 │ Total Volume               88338.2     m³
   9 │ confidence interval lower  82553.6     m³
  10 │ confidence interval upper  94122.7     m³
  11 │ population                     0.945   finite
  12 │ measured plots                11.0     n
  13 │ required plots                 7.0     n
  14 │ missing plots                  0.0     n
  15 │ possible plots               200.0     N
```
"""
function simple_casual_sampling(volume::Vector{<:Real}, plot_area::Real, total_area::Real;
  e::Real=10, α::Real=0.95, lg::Symbol=:en)

  # Validate the language parameter
  if lg ∉ (:pt, :en)
    error("The language 'lg' must be Portuguese ':pt' or English ':en'.")
  end
  # Calculate the total number of possible plots (N)
  N = round(Int, total_area / plot_area)  # N = Total area divided by plot area
  # Number of measured plots (n)
  n = length(volume)
  # Finite population correction factor (f)
  f = 1 - (n / N)  # f = 1 - (sample size / population size)
  # Mean plot volume (x̄)
  x̄ = mean(volume)
  # Coefficient of variation (cv) in percentage
  cv = variation(volume) * 100  # cv = (standard deviation / mean) * 100%
  # Critical t-value for the given confidence level
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)
  # Determine if the population is considered infinite or finite
  if f >= 0.98 || plot_area == 1.0
    # If the correction factor is large or plot area is 1 hectare, consider population infinite
    population = lg == :pt ? "infinita" : "infinite"
    # Variance of the mean (s²x̄) for infinite population
    s²x̄ = var(volume) / n
    # Standard error of the mean (sx̄)
    sx̄ = sqrt(s²x̄)
    # Calculate required number of plots for infinite population
    required_plots = _infinite_population(cv, n, e, α)
  else
    # Population is finite
    population = lg == :pt ? "finita" : "finite"
    # Variance of the mean (s²x̄) with finite population correction
    s²x̄ = (var(volume) / n) * f
    # Standard error of the mean (sx̄)
    sx̄ = sqrt(s²x̄)
    # Calculate required number of plots for finite population
    required_plots = _finite_population(cv, n, N, e, α)
  end
  # Determine the number of additional plots needed
  missing_plots = n > required_plots ? 0 : required_plots - n
  # Generate the inventory report using the calculated statistics
  report = _inventory_report(x̄, cv, Ttab, s²x̄, sx̄, plot_area, f, population, n, required_plots, missing_plots, N, lg)

  return report

end