"""
Calculates the sample size for an infinite population.

# Arguments
- `t::Real`: T-value for the given confidence level.
- `cv::Real`: Coefficient of variation.
- `e::Real`: Desired error margin.

# Returns
- `Float64`: The sample size for an infinite population.
"""
_infinite_sample(t::Real, cv::Real, e::Real) :: Float64 = (t^2 * cv^2) / e^2

"""
Calculates the sample size for a finite population.

# Arguments
- `t::Real`: T-value for the given confidence level.
- `cv::Real`: Coefficient of variation.
- `e::Real`: Desired error margin.
- `N::Int64`: Population size.

# Returns
- `Float64`: The sample size for a finite population.
"""
_finite_sample(t::Real, cv::Real, e::Real, N::Int64) :: Float64 = (t^2 * cv^2) / (e^2 + ((t^2 * cv^2) / N))

"""
Calculates the required sample size for an infinite population using an iterative approach with a maximum iteration limit to prevent infinite loops.

# Arguments
- `cv::Real`: Coefficient of variation.
- `ni::Real`: Initial sample size.
- `e::Real`: Desired error margin.
- `α::Real`: Confidence level.
- `max_iterations::Int=1000`: Maximum number of iterations (default is 1000).

# Returns
- `Int64`: The required sample size for an infinite population.
"""
function _infinite_population(cv::Real, ni::Real, e::Real, α::Real; max_iterations::Int=1000) :: Int64
  iteration = 0
  while iteration < max_iterations
    Ttab1 = -Distributions.quantile(TDist(ni - 1), (1 - α) / 2)
    iter1 = _infinite_sample(Ttab1, cv, e)
    Ttab2 = -Distributions.quantile(TDist(iter1 - 1), (1 - α) / 2)
    iter2 = _infinite_sample(Ttab2, cv, e)
    if round(ni, digits = 3) == round(iter1, digits = 3) || round(ni, digits = 3) == round(iter2, digits = 3)
      return ceil(Integer, ni)
    end
    ni = iter1
    iteration += 1
  end
  error("Maximum number of iterations ($max_iterations) reached without convergence.")
end

"""
Calculates the required sample size for a finite population using an iterative approach with a maximum iteration limit to prevent infinite loops.

# Arguments
- `cv::Real`: Coefficient of variation.
- `ni::Real`: Initial sample size.
- `N::Int`: Population size.
- `e::Real`: Desired error margin.
- `α::Real`: Confidence level.
- `max_iterations::Int=1000`: Maximum number of iterations (default is 1000).

# Returns
- `Int64`: The required sample size for a finite population.
"""
function _finite_population(cv::Real, ni::Real, N::Int, e::Real, α::Real; max_iterations::Int=1000) :: Int64
  iteration = 0
  while iteration < max_iterations
    Ttab1 = -Distributions.quantile(TDist(ni - 1), (1 - α) / 2)
    iter1 = _finite_sample(Ttab1, cv, e, N)
    Ttab2 = -Distributions.quantile(TDist(iter1 - 1), (1 - α) / 2)
    iter2 = _finite_sample(Ttab2, cv, e, N)
    if round(ni, digits = 3) == round(iter1, digits = 3) || round(ni, digits = 3) == round(iter2, digits = 3)
      return ceil(Integer, ni)
    end
    ni = iter1
    iteration += 1
  end
  error("Maximum number of iterations ($max_iterations) reached without convergence.")
end

"""
Performs simple random sampling for forest inventory with specified plot area and total area.

# Arguments
- `volume::Vector{<:Real}`: Volume measurements.
- `plot_area::Real`: Area of each plot.
- `total_area::Real`: Total area.
- `e::Real=10`: Desired error margin (default is 10%).
- `α::Real=0.95`: Confidence level (default is 95%).
- `lg::Symbol=:pt`: Language for the report (default is Portuguese).

# Returns
- `DataFrame`: Inventory report.
"""
function simple_casual_sampling(volume::Vector{<:Real}, plot_area::Real, total_area::Real; e::Real=10, α::Real=0.95, lg::Symbol=:pt)
  if lg ∉ (:pt, :en)
    error("The language 'lg' must be Portuguese 'pt' or English 'en'.")
  end
  N = round(Integer, total_area / plot_area)
  n = length(volume)
  f = 1 - (n / N)
  x̅ = Distributions.mean(volume)
  cv = variation(volume) * 100
  Ttab = -Distributions.quantile(TDist(n - 1), (1 - α) / 2)
  if f >= 0.98 || plot_area == 1.0
    population = lg == :pt ? "infinita" : "infinite"
    s²x̅ = Distributions.var(volume) / n
    sx̅ = √s²x̅
    required_plots = _infinite_population(cv, n, e, α)
  else
    population = lg == :pt ? "finita" : "finite"
    s²x̅ = (Distributions.var(volume) / n) * f
    sx̅ = √s²x̅
    required_plots = _finite_population(cv, n, N, e, α)
  end
  missing_plots = n > required_plots ? 0 : required_plots - n
  _inventory_report(x̅, cv, Ttab, s²x̅, sx̅, plot_area, f, population, n, required_plots, missing_plots, N, lg)
end

"""
Performs simple random sampling for forest inventory with total area.

# Arguments
- `volume::Vector{<:Real}`: Volume measurements.
- `total_area::Real`: Total area.
- `e::Real=10`: Desired error margin (default is 10%).
- `α::Real=0.95`: Confidence level (default is 95%).
- `lg::Symbol=:pt`: Language for the report (":pt" for Portuguese - default, ":en" for English).

# Returns
- `DataFrame`: Inventory report.
"""
function simple_casual_sampling(volume::Vector{<:Real}, total_area::Real; e::Real=10, α::Real=0.95, lg::Symbol=:pt)
  if lg ∉ (:pt, :en)
    error("The language 'lg' must be Portuguese 'pt' ou English 'en'.")
  end
  N = total_area
  n = length(volume)
  x̅ = Distributions.mean(volume)
  cv = variation(volume) * 100
  Ttab = -Distributions.quantile(TDist(n - 1), (1 - α) / 2)
  f = NaN64
  population = lg == :pt ? "infinita" : "infinite"
  s²x̅ = Distributions.var(volume) / n
  sx̅ = √s²x̅
  required_plots = _infinite_population(cv, n, e, α)
  missing_plots = n > required_plots ? 0 : required_plots - n
  _inventory_report(x̅, cv, Ttab, s²x̅, sx̅, "u.a.", f, population, n, required_plots, missing_plots, N, lg)
end