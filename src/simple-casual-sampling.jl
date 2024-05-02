_infinite_sample(t::Real, cv::Real, e::Real) :: Float64 = (t^2 * cv^2) / e^2

_finite_sample(t::Real, cv::Real, e::Real, N::Int64) :: Float64 = (t^2 * cv^2) / (e^2 + ((t^2 * cv^2) / N))

function _infinite_population(cv::Real, ni::Real, e::Real, α::Real) :: Int64
  Ttab1 = -quantile(TDist(ni - 1), (1 - α) / 2)
  iter1 = _infinite_sample(Ttab1, cv, e)
  Ttab2 = -quantile(TDist(iter1 - 1), (1 - α) / 2)
  iter2 = _infinite_sample(Ttab2, cv, e)
  if round(ni, digits = 3) != round(iter1, digits = 3) && round(ni, digits = 3) != round(iter2, digits = 3)
    _infinite_population(cv, iter1, e, α)
  else
    return ceil(Integer, ni)
  end
end

function _finite_population(cv::Real, ni::Real, N::Int, e::Real, α::Real) :: Int64
  Ttab1 = -quantile(TDist(ni - 1), (1 - α) / 2)
  iter1 = _finite_sample(Ttab1, cv, e, N)
  Ttab2 = -quantile(TDist(iter1 - 1), (1 - α) / 2)
  iter2 = _finite_sample(Ttab2, cv, e, N)
  if round(ni, digits = 3) != round(iter1, digits = 3) && round(ni, digits = 3) != round(iter2, digits = 3)
    _finite_population(cv, iter1, N, e, α)
  else
    return ceil(Integer, ni)
  end
end

function simple_casual_sampling(volume::Vector{<:Real}, plot_area::Real, total_area::Real; e::Real=10, α::Real=0.95, lg::Symbol=:pt)
  if lg ∉ (:pt, :en)
    error("The language 'lg' must be Portuguese 'pt' or English 'en'.")
  end
  N = round(Integer, total_area / plot_area)
  n = length(volume)
  f = 1 - (n /  N)
  x̅ = mean(volume)
  cv = variation(volume) * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)
  if f >= 0.98 || plot_area == 1.0
    population = lg == :pt ? "infinita" : "infinite"
    s²x̅ = var(volume) / n
    sx̅ = √s²x̅
    required_plots = _infinite_population(cv, n, e, α)
  else
    population = lg == :pt ? "finita" : "finite"
    s²x̅ = (var(volume) / n) * f
    sx̅ = √s²x̅
    required_plots = _finite_population(cv, n, N, e, α)
  end
  missing_plots = n > required_plots ? 0 : required_plots - n
  _inventory_report(x̅, cv, Ttab, s²x̅, sx̅, plot_area, f, population, n, required_plots, missing_plots, N, lg)
end

function simple_casual_sampling(volume::Vector{<:Real}, total_area::Real; e::Real=10, α::Real=0.95, lg::Symbol=:pt)
  if lg ∉ (:pt, :en)
    error("The language 'lg' must be Portuguese 'pt' or English 'en'.")
  end
  N = total_area
  n = length(volume)
  x̅ = mean(volume)
  cv = variation(volume) * 100
  Ttab = -quantile(TDist(n - 1), (1 - α) / 2)
  f = NaN64
  population = lg == :pt ? "infinita" : "infinite"
  s²x̅ = var(volume) / n
  sx̅ = √s²x̅
  required_plots = _infinite_population(cv, n, e, α)
  missing_plots = n > required_plots ? 0 : required_plots - n
  _inventory_report(x̅, cv, Ttab, s²x̅, sx̅, "u.a.", f, population, n, required_plots, missing_plots, N, lg)
end