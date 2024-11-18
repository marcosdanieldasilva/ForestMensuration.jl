# Calculates the number of classes using Sturges' formula.
# Arguments
# - `n::Int`: The number of observations.
# Returns
# - `Int`: The number of classes.
_sturges(n::Int) = ceil(Int, log2(n)) + 1


# Calculates the class center for a given value and class width.
# Arguments
# - `x::Real`: The value.
# - `hi::Real`: The class width.
# Returns
# - `Real`: The class center.
_class_center(x::Real, hi::Real) = round(x / hi) * hi + (hi / 2)


# Calculates the amplitude (range) of a vector of values.
# Arguments
# - `x::Vector`: The vector of values.
# Returns
# - `Real`: The amplitude (range).
_amplitude(x::Vector) = maximum(x) - minimum(x)


# Calculates the class breadth (width) for a given amplitude and number of classes.
# Arguments
# - `h::Real`: The amplitude.
# - `k::Int`: The number of classes.
# Returns
# - `Real`: The class breadth (width).
function _class_breadth(h::Real, k::Int)
  hi = h / k
  log_hi = log10(hi)
  if log_hi >= 0
    step = exp10(floor(log_hi))
    rate = hi / step
    if rate <= 1.1
      nothing
    elseif rate <= 2.2
      step *= 2
    elseif rate <= 5.5
      step *= 5
    else
      step *= 10
    end
    return step
  else
    step = exp10(-floor(log_hi))
    rate = hi * step
    if rate <= 1.1
      nothing
    elseif rate <= 2.2
      step /= 2
    elseif rate <= 5.5
      step /= 5
    else
      step /= 10
    end
    return step^-1
  end
end


# Calculates the simple frequency of unique values in a vector.
# Arguments
# - `x::Vector`: The vector of values.
# Returns
# - `Vector`: The frequency of each unique value.
_simple_frequency(x::Vector) = map(i -> count(==(i), x), unique(x) |> sort)

"""
  frequency_table(x::Vector{<:Real}, hi::Real)

Creates a frequency table for a vector of values given a class width.

# Arguments
- `x::Vector{<:Real}`: The vector of values.
- `hi::Real`: The class width.

# Returns
- `DataFrame`: A DataFrame containing the frequency table.
"""
function frequency_table(x::Vector{<:Real}, hi::Real)
  n = length(x)
  cc = _class_center.(x, hi)
  Xi = unique(cc) |> sort
  LI = Xi .- hi / 2
  LS = Xi .+ hi / 2
  fi = _simple_frequency(cc)
  Fi = cumsum(fi)
  fri = (fi ./ n) .* 100
  Fri = cumsum(fri)
  DataFrame(LI=LI, Xi=Xi, LS=LS, fi=fi, Fi=Fi, fri=fri, Fri=Fri)
end

"""
  frequency_table(x::Vector{<:Real})

Creates a frequency table for a vector of values.

# Arguments
- `x::Vector{<:Real}`: The vector of values.

# Returns
- `DataFrame`: A DataFrame containing the frequency table.
"""
function frequency_table(x::Vector{<:Real})
  h = _amplitude(x)
  k = _sturges(length(x))
  hi = _class_breadth(h, k)
  frequency_table(x, hi)
end

"""
  frequency_table(g::Symbol, x::Symbol, data::AbstractDataFrame)

Creates a frequency table for grouped data in a DataFrame.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `x::S`: The symbol representing the values variable.
- `data::AbstractDataFrame`: The DataFrame containing the data.

# Returns
- `DataFrame`: A DataFrame containing the frequency table for each group.
"""
function frequency_table(g::Symbol, x::Symbol, data::AbstractDataFrame)
  combine(groupby(data, g)) do df
    frequency_table(df[:, x])
  end
end

"""
  frequency_table(g::Symbol, x::Symbol, hi::Real, data::AbstractDataFrame)

Creates a frequency table for grouped data in a DataFrame with a specified class width.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `x::S`: The symbol representing the values variable.
- `hi::Real`: The class width.
- `data::AbstractDataFrame`: The DataFrame containing the data.

# Returns
- `DataFrame`: A DataFrame containing the frequency table for each group.
"""
function frequency_table(g::Symbol, x::Symbol, hi::Real, data::AbstractDataFrame)
  combine(groupby(data, g)) do df
    frequency_table(df[:, x], hi)
  end
end

"""
  diametric_table(x::Vector{<:Real}, hi::Real; plot_area::Real=1.0)

Creates a diametric table for a vector of values given a class width and plot area.

# Arguments
- `x::Vector{<:Real}`: The vector of values.
- `hi::Real`: The class width.
- `plot_area::Real=1.0`: The plot area in hectares (default is 1.0).

# Returns
- `DataFrame`: A DataFrame containing the diametric table.
"""
function diametric_table(x::Vector{<:Real}, hi::Real; plot_area::Real=1.0)
  ftable = frequency_table(x, hi)
  ftable.g = (π / 40000) .* ftable.Xi .^ 2
  ftable.ng = ftable.fi .* ftable.g
  ftable.∑ng = cumsum(ftable.ng)
  if plot_area != 1.0
    ftable.fi_ha = ftable.fi * (1 / plot_area)
    ftable.Fi_ha = cumsum(ftable.fi_ha)
    ftable.ng_ha = ftable.ng * (1 / plot_area)
    ftable.∑ng_ha = cumsum(ftable.ng_ha)
  end
  return ftable
end

"""
  diametric_table(x::Vector{<:Real}; plot_area::Real=1.0

Creates a diametric table for a vector of values.

# Arguments
- `x::Vector{<:Real}`: The vector of values.
- `plot_area::Real=1.0`: The plot area in hectares (default is 1.0).

# Returns
- `DataFrame`: A DataFrame containing the diametric table.
"""
function diametric_table(x::Vector{<:Real}; plot_area::Real=1.0)
  h = _amplitude(x)
  k = _sturges(length(x))
  hi = _class_breadth(h, k)
  diametric_table(x, hi, plot_area=plot_area)
end

"""
  diametric_table(g::Symbol, x::Symbol, data::AbstractDataFrame; plot_area::Real=1.0)

Creates a diametric table for grouped data in a DataFrame.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `x::S`: The symbol representing the values variable.
- `data::AbstractDataFrame`: The DataFrame containing the data.
- `plot_area::Real=1.0`: The plot area in hectares (default is 1.0).

# Returns
- `DataFrame`: A DataFrame containing the diametric table for each group.
"""
function diametric_table(g::Symbol, x::Symbol, data::AbstractDataFrame; plot_area::Real=1.0)
  combine(groupby(data, g)) do df
    diametric_table(df[:, x], plot_area=plot_area)
  end
end

"""
  diametric_table(g::Symbol, x::Symbol, hi::Real, data::AbstractDataFrame; plot_area::Real=1.0)
  
Creates a diametric table for grouped data in a DataFrame with a specified class width.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `x::S`: The symbol representing the values variable.
- `hi::Real`: The class width.
- `data::AbstractDataFrame`: The DataFrame containing the data.
- `plot_area::Real=1.0`: The plot area in hectares (default is 1.0).

# Returns
- `DataFrame`: A DataFrame containing the diametric table for each group.
"""
function diametric_table(g::Symbol, x::Symbol, hi::Real, data::AbstractDataFrame; plot_area::Real=1.0)
  combine(groupby(data, g)) do df
    diametric_table(df[:, x], hi, plot_area=plot_area)
  end
end