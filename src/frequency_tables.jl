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

# Example
```julia

# Define the vector of values
julia> x = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];

# Calculate the frequency table with The class width of 2
julia> frequency_table(x, 2)
8×7 DataFrame
 Row │ LI       Xi       LS       fi     Fi     fri      Fri     
     │ Float64  Float64  Float64  Int64  Int64  Float64  Float64
─────┼───────────────────────────────────────────────────────────
   1 │    10.0     11.0     12.0      1      1     10.0     10.0
   2 │    12.0     13.0     14.0      1      2     10.0     20.0
   3 │    14.0     15.0     16.0      1      3     10.0     30.0
   4 │    16.0     17.0     18.0      2      5     20.0     50.0
   5 │    18.0     19.0     20.0      1      6     10.0     60.0
   6 │    20.0     21.0     22.0      2      8     20.0     80.0
   7 │    22.0     23.0     24.0      1      9     10.0     90.0
   8 │    24.0     25.0     26.0      1     10     10.0    100.0
```
"""
function frequency_table(x::Vector{<:Real}, hi::Real)
  if hi <= 0
    throw(DomainError("The class width must be positive."))
  end
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

# Example
```julia

# Define the vector of values
julia> x = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];

# Calculate the frequency table with auto class width
julia> frequency_table(x)
4×7 DataFrame
 Row │ LI       Xi       LS       fi     Fi     fri      Fri     
     │ Float64  Float64  Float64  Int64  Int64  Float64  Float64
─────┼───────────────────────────────────────────────────────────
   1 │    10.0     12.5     15.0      2      2     20.0     20.0
   2 │    15.0     17.5     20.0      3      5     30.0     50.0
   3 │    20.0     22.5     25.0      4      9     40.0     90.0
   4 │    25.0     27.5     30.0      1     10     10.0    100.0
```
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

# Example
```julia
julia> using DataFrames

# Define the data'
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];
julia> species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"];
julia> data = DataFrame(species=species, diameters=diameters);

# Calculate the frequency table by group with auto class width
julia> frequency_table(:species, :diameters, data)
8×8 DataFrame
 Row │ species  LI       Xi       LS       fi     Fi     fri      Fri     
     │ String   Float64  Float64  Float64  Int64  Int64  Float64  Float64
─────┼────────────────────────────────────────────────────────────────────
   1 │ Oak         10.0     11.0     12.0      1      1     20.0     20.0
   2 │ Oak         12.0     13.0     14.0      1      2     20.0     40.0
   3 │ Oak         14.0     15.0     16.0      1      3     20.0     60.0
   4 │ Oak         16.0     17.0     18.0      2      5     40.0    100.0
   5 │ Pine        18.0     19.0     20.0      1      1     20.0     20.0
   6 │ Pine        20.0     21.0     22.0      2      3     40.0     60.0
   7 │ Pine        22.0     23.0     24.0      1      4     20.0     80.0
   8 │ Pine        24.0     25.0     26.0      1      5     20.0    100.0
```
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

# Example
```julia
julia> using DataFrames

# Define the data'
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];
julia> species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"];
julia> data = DataFrame(species=species, diameters=diameters);

# Calculate the frequency table by group with with The class width of 2
julia> frequency_table(:species, :diameters, 2, data)
8×8 DataFrame
 Row │ species  LI       Xi       LS       fi     Fi     fri      Fri     
     │ String   Float64  Float64  Float64  Int64  Int64  Float64  Float64
─────┼────────────────────────────────────────────────────────────────────
   1 │ Oak         10.0     11.0     12.0      1      1     20.0     20.0
   2 │ Oak         12.0     13.0     14.0      1      2     20.0     40.0
   3 │ Oak         14.0     15.0     16.0      1      3     20.0     60.0
   4 │ Oak         16.0     17.0     18.0      2      5     40.0    100.0
   5 │ Pine        18.0     19.0     20.0      1      1     20.0     20.0
   6 │ Pine        20.0     21.0     22.0      2      3     40.0     60.0
   7 │ Pine        22.0     23.0     24.0      1      4     20.0     80.0
   8 │ Pine        24.0     25.0     26.0      1      5     20.0    100.0
```
"""
function frequency_table(g::Symbol, x::Symbol, hi::Real, data::AbstractDataFrame)
  combine(groupby(data, g)) do df
    frequency_table(df[:, x], hi)
  end
end

"""
    diametric_table(d::Vector{<:Real}, hi::Real; plot_area::Real=1.0)

Creates a diametric table for a vector of diameter values given a class width and plot area.

# Arguments
- `d::Vector{<:Real}`: The vector of diameter values.
- `hi::Real`: The class width.
- `plot_area::Real=1.0`: The plot area in hectares (default is 1.0).

# Returns
- `DataFrame`: A DataFrame containing the diametric table.

# Example
```julia

# Define the vector of diameter values
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];

# Calculate the diametric table with with The class width of 2 and plot area of 0.05 ha (500 m²)
julia> diametric_table(diameters, 2, plot_area=0.05)
8×14 DataFrame
 Row │ LI       Xi       LS       fi     Fi     fri      Fri      g           ng          ∑ng         fi_ha    Fi_ha    ng_ha     ∑ng_ha   
     │ Float64  Float64  Float64  Int64  Int64  Float64  Float64  Float64     Float64     Float64     Float64  Float64  Float64   Float64  
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │    10.0     11.0     12.0      1      1     10.0     10.0  0.00950332  0.00950332  0.00950332     20.0     20.0  0.190066  0.190066
   2 │    12.0     13.0     14.0      1      2     10.0     20.0  0.0132732   0.0132732   0.0227765      20.0     40.0  0.265465  0.455531
   3 │    14.0     15.0     16.0      1      3     10.0     30.0  0.0176715   0.0176715   0.040448       20.0     60.0  0.353429  0.80896
   4 │    16.0     17.0     18.0      2      5     20.0     50.0  0.022698    0.045396    0.085844       40.0    100.0  0.90792   1.71688
   5 │    18.0     19.0     20.0      1      6     10.0     60.0  0.0283529   0.0283529   0.114197       20.0    120.0  0.567057  2.28394
   6 │    20.0     21.0     22.0      2      8     20.0     80.0  0.0346361   0.0692721   0.183469       40.0    160.0  1.38544   3.66938
   7 │    22.0     23.0     24.0      1      9     10.0     90.0  0.0415476   0.0415476   0.225017       20.0    180.0  0.830951  4.50033
   8 │    24.0     25.0     26.0      1     10     10.0    100.0  0.0490874   0.0490874   0.274104       20.0    200.0  0.981748  5.48208
```
"""
function diametric_table(d::Vector{<:Real}, hi::Real; plot_area::Real=1.0)
  if hi <= 0
    throw(DomainError("The class width must be positive."))
  elseif plot_area <= 0
    throw(DomainError("The plot area must be positive."))
  elseif any(x -> x <= 0, d)
    throw(DomainError("Diameters must be positive"))
  else
    ftable = frequency_table(d, hi)
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
end

"""
    diametric_table(d::Vector{<:Real}; plot_area::Real=1.0)

Creates a diametric table for a vector of diameter values.

# Arguments
- `d::Vector{<:Real}`: The vector of diameter values.
- `plot_area::Real=1.0`: The plot area in hectares (default is 1.0).

# Returns
- `DataFrame`: A DataFrame containing the diametric table.

# Example
```julia

# Define the vector of diameter values
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];

# Calculate the diametric table with with auto class width and plot area of 0.05 ha (500 m²)
4×14 DataFrame
 Row │ LI       Xi       LS       fi     Fi     fri      Fri      g          ng         ∑ng        fi_ha    Fi_ha    ng_ha     ∑ng_ha   
     │ Float64  Float64  Float64  Int64  Int64  Float64  Float64  Float64    Float64    Float64    Float64  Float64  Float64   Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │    10.0     12.5     15.0      2      2     20.0     20.0  0.0122718  0.0245437  0.0245437     40.0     40.0  0.490874  0.490874
   2 │    15.0     17.5     20.0      3      5     30.0     50.0  0.0240528  0.0721585  0.0967021     60.0    100.0  1.44317   1.93404
   3 │    20.0     22.5     25.0      4      9     40.0     90.0  0.0397608  0.159043   0.255745      80.0    180.0  3.18086   5.11491
   4 │    25.0     27.5     30.0      1     10     10.0    100.0  0.0593957  0.0593957  0.315141      20.0    200.0  1.18791   6.30282
```
"""
function diametric_table(d::Vector{<:Real}; plot_area::Real=1.0)
  h = _amplitude(d)
  k = _sturges(length(d))
  hi = _class_breadth(h, k)
  diametric_table(d, hi, plot_area=plot_area)
end

"""
    diametric_table(g::Symbol, d::Symbol, data::AbstractDataFrame; plot_area::Real=1.0)

Creates a diametric table for grouped data in a DataFrame.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `d::S`: The symbol representing the diameter values variable.
- `data::AbstractDataFrame`: The DataFrame containing the data.
- `plot_area::Real=1.0`: The plot area in hectares (default is 1.0).

# Returns
- `DataFrame`: A DataFrame containing the diametric table for each group.

# Example
```julia
julia> using DataFrames

# Define the data'
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];
julia> species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"];
julia> data = DataFrame(species=species, diameters=diameters);

# Calculate the diametric table by group with with auto class width and plot area of 0.05 ha (500 m²)
julia> diametric_table(:species, :diameters, data, plot_area=0.05)
8×15 DataFrame
 Row │ species  LI       Xi       LS       fi     Fi     fri      Fri      g           ng          ∑ng         fi_ha    Fi_ha    ng_ha     ∑ng_ha   
     │ String   Float64  Float64  Float64  Int64  Int64  Float64  Float64  Float64     Float64     Float64     Float64  Float64  Float64   Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Oak         10.0     11.0     12.0      1      1     20.0     20.0  0.00950332  0.00950332  0.00950332     20.0     20.0  0.190066  0.190066
   2 │ Oak         12.0     13.0     14.0      1      2     20.0     40.0  0.0132732   0.0132732   0.0227765      20.0     40.0  0.265465  0.455531
   3 │ Oak         14.0     15.0     16.0      1      3     20.0     60.0  0.0176715   0.0176715   0.040448       20.0     60.0  0.353429  0.80896
   4 │ Oak         16.0     17.0     18.0      2      5     40.0    100.0  0.022698    0.045396    0.085844       40.0    100.0  0.90792   1.71688
   5 │ Pine        18.0     19.0     20.0      1      1     20.0     20.0  0.0283529   0.0283529   0.0283529      20.0     20.0  0.567057  0.567057
   6 │ Pine        20.0     21.0     22.0      2      3     40.0     60.0  0.0346361   0.0692721   0.097625       40.0     60.0  1.38544   1.9525
   7 │ Pine        22.0     23.0     24.0      1      4     20.0     80.0  0.0415476   0.0415476   0.139173       20.0     80.0  0.830951  2.78345
   8 │ Pine        24.0     25.0     26.0      1      5     20.0    100.0  0.0490874   0.0490874   0.18826        20.0    100.0  0.981748  3.7652
```
"""
function diametric_table(g::Symbol, d::Symbol, data::AbstractDataFrame; plot_area::Real=1.0)
  combine(groupby(data, g)) do df
    diametric_table(df[:, d], plot_area=plot_area)
  end
end

"""
    diametric_table(g::Symbol, d::Symbol, hi::Real, data::AbstractDataFrame; plot_area::Real=1.0)
  
Creates a diametric table for grouped data in a DataFrame with a specified class width.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `d::S`: The symbol representing the diameter values variable.
- `hi::Real`: The class width.
- `data::AbstractDataFrame`: The DataFrame containing the data.
- `plot_area::Real=1.0`: The plot area in hectares (default is 1.0).

# Returns
- `DataFrame`: A DataFrame containing the diametric table for each group.

# Example
```julia
julia> using DataFrames

# Define the data'
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];
julia> species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"];
julia> data = DataFrame(species=species, diameters=diameters);

# Calculate the diametric table by group with with The class width of 3 and plot area of 0.05 ha (500 m²)
julia> diametric_table(:species, :diameters, 3, data, plot_area=0.05)
6×15 DataFrame
 Row │ species  LI       Xi       LS       fi     Fi     fri      Fri      g          ng         ∑ng        fi_ha    Fi_ha    ng_ha     ∑ng_ha   
     │ String   Float64  Float64  Float64  Int64  Int64  Float64  Float64  Float64    Float64    Float64    Float64  Float64  Float64   Float64
─────┼───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ Oak         12.0     13.5     15.0      3      3     60.0     60.0  0.0143139  0.0429416  0.0429416     60.0     60.0  0.858833  0.858833
   2 │ Oak         15.0     16.5     18.0      1      4     20.0     80.0  0.0213825  0.0213825  0.0643241     20.0     80.0  0.427649  1.28648
   3 │ Oak         18.0     19.5     21.0      1      5     20.0    100.0  0.0298648  0.0298648  0.0941889     20.0    100.0  0.597295  1.88378
   4 │ Pine        18.0     19.5     21.0      2      2     40.0     40.0  0.0298648  0.0597295  0.0597295     40.0     40.0  1.19459   1.19459
   5 │ Pine        21.0     22.5     24.0      1      3     20.0     60.0  0.0397608  0.0397608  0.0994903     20.0     60.0  0.795216  1.98981
   6 │ Pine        24.0     25.5     27.0      2      5     40.0    100.0  0.0510705  0.102141   0.201631      40.0    100.0  2.04282   4.03263
```
"""
function diametric_table(g::Symbol, d::Symbol, hi::Real, data::AbstractDataFrame; plot_area::Real=1.0)
  combine(groupby(data, g)) do df
    diametric_table(df[:, d], hi, plot_area=plot_area)
  end
end