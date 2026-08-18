# Calculates the number of classes using Sturges' formula.
_sturges(n::Int) = ceil(Int, log2(n)) + 1
# Calculates the class center for a given value and class width.
_classcenter(x::Real, hi::Real) = round(x / hi) * hi + (hi / 2)
# Calculates the amplitude (range) of a vector of values.
_amplitude(x::AbstractVector) = maximum(x) - minimum(x)
# Calculates the class breadth (width) for a given amplitude and number of classes.
function _classbreadth(h::Real, k::Int)
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
_simplefrequency(x::AbstractVector) = map(i -> count(==(i), x), unique(x) |> sort)

# Unit in which a bare (unitless) class width is interpreted for a unitful sample:
# the same unit the observations themselves are expressed in.
_sampleunit(x::AbstractVector{<:Quantity}) = unit(first(x))

"""
    frequencytable(x::AbstractVector{<:Real}, hi::Real)
    frequencytable(x::AbstractVector{<:Quantity}, hi::Union{Quantity,Real})

Creates a frequency table for a vector of values given a class width.

# Arguments
- `x::AbstractVector{<:Real}`: The vector of values.
- `hi::Real`: The class width.

# Units

`frequencytable` is unit-agnostic: it classifies whatever variable it is given
(diameters, heights, volumes, ages). When `x` carries `Unitful` quantities, the class
limits `LI`, `Xi` and `LS` are returned in the same unit as `x`, and a class width given
as a plain number is interpreted in that unit as well. When `x` holds plain numbers the
table is returned unitless — no unit is assumed, because the variable is unknown. Use
[`diametrictable`](@ref) for diameter-specific tables, where a plain number is taken to
be centimeters.

# Returns
- `DataFrame`: A DataFrame containing the frequency table.

# Example
```julia-repl

# Define the vector of values
julia> x = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];

# Calculate the frequency table with The class width of 2
julia> frequencytable(x, 2)
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

# the same sample carrying units keeps them in the class limits
julia> frequencytable([10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]u"cm", 2u"cm")
8×7 DataFrame
 Row │ LI        Xi        LS        fi     Fi     fri      Fri
     │ Quantity  Quantity  Quantity  Int64  Int64  Float64  Float64
─────┼────────────────────────────────────────────────────────────────
   1 │  10.0 cm   11.0 cm   12.0 cm      1      1     10.0     10.0
  ⋮  │    ⋮         ⋮         ⋮        ⋮      ⋮       ⋮        ⋮
```
"""
function frequencytable(x::AbstractVector{<:Real}, hi::Real)
  if hi <= 0
    throw(DomainError("The class width must be positive."))
  end
  n = length(x)
  cc = _classcenter.(x, hi)
  Xi = unique(cc) |> sort
  LI = Xi .- hi / 2
  LS = Xi .+ hi / 2
  fi = _simplefrequency(cc)
  Fi = cumsum(fi)
  fri = (fi ./ n) .* 100
  Fri = cumsum(fri)
  DataFrame(LI=LI, Xi=Xi, LS=LS, fi=fi, Fi=Fi, fri=fri, Fri=Fri)
end

function frequencytable(x::AbstractVector{<:Quantity}, hi::Union{Quantity,Real})
  ux = _sampleunit(x)
  hiq = withunit(hi, ux)
  # classification happens on plain numbers in a single unit, then the class limits are
  # restored so the table speaks the same unit as the observations
  ftable = frequencytable(ustrip.(uconvert.(ux, x)), ustrip(uconvert(ux, hiq)))
  ftable.LI = ftable.LI * ux
  ftable.Xi = ftable.Xi * ux
  ftable.LS = ftable.LS * ux
  return ftable
end

"""
    frequencytable(x::AbstractVector{<:Real})
    frequencytable(x::AbstractVector{<:Quantity})

Creates a frequency table for a vector of values, choosing the class width automatically
from Sturges' rule.

# Arguments
- `x::AbstractVector{<:Real}`: The vector of values.

# Units

The class limits are returned in the unit of `x` when it carries `Unitful` quantities,
and unitless otherwise.

# Returns
- `DataFrame`: A DataFrame containing the frequency table.

# Example
```julia-repl

# Define the vector of values
julia> x = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];

# Calculate the frequency table with auto class width
julia> frequencytable(x)
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
function frequencytable(x::AbstractVector{<:Real})
  h = _amplitude(x)
  k = _sturges(length(x))
  hi = _classbreadth(h, k)
  frequencytable(x, hi)
end

function frequencytable(x::AbstractVector{<:Quantity})
  ux = _sampleunit(x)
  h = ustrip(uconvert(ux, _amplitude(x)))
  k = _sturges(length(x))
  hi = _classbreadth(h, k)
  frequencytable(x, hi * ux)
end

"""
    frequencytable(g::Symbol, x::Symbol, data::AbstractDataFrame)

Creates a frequency table for grouped data in a DataFrame.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `x::S`: The symbol representing the values variable.
- `data::AbstractDataFrame`: The DataFrame containing the data.

# Units

The value column may hold either plain numbers or `Unitful` quantities; the class limits
are returned in the same form.

# Returns
- `DataFrame`: A DataFrame containing the frequency table for each group.

# Example
```julia-repl
julia> using DataFrames

# Define the data'
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];
julia> species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"];
julia> data = DataFrame(species=species, diameters=diameters);

# Calculate the frequency table by group with auto class width
julia> frequencytable(:species, :diameters, data)
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
function frequencytable(g::Symbol, x::Symbol, data::AbstractDataFrame)
  combine(groupby(data, g)) do df
    frequencytable(df[:, x])
  end
end

"""
    frequencytable(g::Symbol, x::Symbol, hi::Union{Quantity,Real}, data::AbstractDataFrame)

Creates a frequency table for grouped data in a DataFrame with a specified class width.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `x::S`: The symbol representing the values variable.
- `hi::Real`: The class width. When the value column carries units, a plain number is interpreted in the unit of that column.
- `data::AbstractDataFrame`: The DataFrame containing the data.

# Returns
- `DataFrame`: A DataFrame containing the frequency table for each group.

# Example
```julia-repl
julia> using DataFrames

# Define the data'
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];
julia> species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"];
julia> data = DataFrame(species=species, diameters=diameters);

# Calculate the frequency table by group with with The class width of 2
julia> frequencytable(:species, :diameters, 2, data)
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
function frequencytable(g::Symbol, x::Symbol, hi::Union{Quantity,Real}, data::AbstractDataFrame)
  combine(groupby(data, g)) do df
    frequencytable(df[:, x], hi)
  end
end

"""
    diametrictable(d::AbstractVector{<:Len}, hi::Union{Len,Real}; plot_area::Union{Area,Real}=1.0u"ha")
    diametrictable(d::AbstractVector{<:Real}, hi::Real; plot_area::Real=1.0)

Creates a diametric table for a vector of diameter values given a class width and plot area.

# Arguments
- `d::AbstractVector{<:Len}`: The vector of diameter values.
- `hi::Len`: The class width.
- `plot_area::Area=1.0u"ha"`: The sampled plot area.

# Units

Unlike the general-purpose [`frequencytable`](@ref), this function is diameter-specific,
so plain numbers are promoted with the package defaults: diameters and class widths in
centimeters (`cm`) and the plot area in hectares (`ha`). The returned table always carries
units:

- `LI`, `Xi`, `LS` in the unit of `d`;
- `g`, `ng` and `∑ng` in square feet (`ft^2`) when `d` is imperial and square meters (`m^2`) otherwise;
- `fi_ha` and `Fi_ha` per acre (`ac^-1`) when `d` is imperial and per hectare (`ha^-1`) otherwise, with `ng_ha` and `∑ng_ha` following the same reference area.

The per-area columns are omitted when the plot already covers exactly one reference area
(1 ha in the metric system, 1 ac in the imperial one), since the expansion would be a
no-op.

# Returns
- `DataFrame`: A DataFrame containing the diametric table.

# Example
```julia-repl

# Define the vector of diameter values
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]u"cm";

# Calculate the diametric table with The class width of 2 cm and plot area of 500 m² (0.05 ha)
julia> diametrictable(diameters, 2u"cm", plot_area=500u"m^2")
8×14 DataFrame
 Row │ LI       Xi       LS       fi     Fi     fri      Fri      g              ng             ∑ng            fi_ha       Fi_ha        ng_ha              ∑ng_ha
     │ Quantity Quantity Quantity Int64  Int64  Float64  Float64  Quantity       Quantity       Quantity       Quantity    Quantity     Quantity           Quantity
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 10.0 cm  11.0 cm  12.0 cm      1      1     10.0     10.0  0.0095 m^2     0.0095 m^2     0.0095 m^2     20.0 ha^-1  20.0 ha^-1   0.190066 m^2 ha^-1  0.190066 m^2 ha^-1
  ⋮  │   ⋮        ⋮        ⋮        ⋮      ⋮       ⋮        ⋮          ⋮              ⋮              ⋮             ⋮           ⋮              ⋮                  ⋮

# the same call without units: diameters and class width in cm, plot area in ha
julia> diametrictable([10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0], 2, plot_area=0.05)
```
"""
function diametrictable(d::AbstractVector{<:Len}, hi::Union{Len,Real}; plot_area::Union{Area,Real}=1.0u"ha")
  ud = _sampleunit(d)
  hi = withunit(hi, ud)
  plot_area = withunit(plot_area, AUNIT)
  if hi <= zero(hi)
    throw(DomainError("The class width must be positive."))
  elseif plot_area <= zero(plot_area)
    throw(DomainError("The plot area must be positive."))
  elseif any(x -> x <= zero(x), d)
    throw(DomainError("Diameters must be positive"))
  end
  ftable = frequencytable(d, hi)
  ftable.g = basalarea.(ftable.Xi)
  ftable.ng = ftable.fi .* ftable.g
  ftable.∑ng = cumsum(ftable.ng)
  # expansion factor to the reference area of the system the diameters were measured in
  EF = expansionfactor(ud, plot_area)
  if !isone(ustrip(EF))
    ftable.fi_ha = ftable.fi * EF
    ftable.Fi_ha = cumsum(ftable.fi_ha)
    ftable.ng_ha = ftable.ng * EF
    ftable.∑ng_ha = cumsum(ftable.ng_ha)
  end
  return ftable
end

diametrictable(d::AbstractVector{<:Real}, hi::Union{Len,Real}; plot_area::Union{Area,Real}=1.0) =
  diametrictable(d * DUNIT, hi; plot_area=withunit(plot_area, AUNIT))

"""
    diametrictable(d::AbstractVector{<:Len}; plot_area::Union{Area,Real}=1.0u"ha")
    diametrictable(d::AbstractVector{<:Real}; plot_area::Real=1.0)

Creates a diametric table for a vector of diameter values, choosing the class width
automatically from Sturges' rule.

# Arguments
- `d::AbstractVector{<:Len}`: The vector of diameter values.
- `plot_area::Area=1.0u"ha"`: The sampled plot area.

# Units

Plain numbers are taken to be diameters in centimeters (`cm`) and a plot area in hectares
(`ha`). See the method with an explicit class width for the units of every output column.

# Returns
- `DataFrame`: A DataFrame containing the diametric table.

# Example
```julia-repl

# Define the vector of diameter values
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]u"cm";

# Calculate the diametric table with auto class width and plot area of 500 m² (0.05 ha)
julia> diametrictable(diameters, plot_area=500u"m^2")
4×14 DataFrame
 Row │ LI       Xi       LS       fi     Fi     fri      Fri      g            ng           ∑ng          fi_ha       Fi_ha       ng_ha                ∑ng_ha
     │ Quantity Quantity Quantity Int64  Int64  Float64  Float64  Quantity     Quantity     Quantity     Quantity    Quantity    Quantity             Quantity
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 10.0 cm  12.5 cm  15.0 cm      2      2     20.0     20.0  0.0122718 m^2 0.0245437 m^2 0.0245437 m^2 40.0 ha^-1 40.0 ha^-1 0.490874 m^2 ha^-1  0.490874 m^2 ha^-1
  ⋮  │   ⋮        ⋮        ⋮        ⋮      ⋮       ⋮        ⋮          ⋮            ⋮            ⋮            ⋮          ⋮               ⋮                   ⋮
```
"""
function diametrictable(d::AbstractVector{<:Len}; plot_area::Union{Area,Real}=1.0u"ha")
  ud = _sampleunit(d)
  h = ustrip(uconvert(ud, _amplitude(d)))
  k = _sturges(length(d))
  hi = _classbreadth(h, k)
  diametrictable(d, hi * ud, plot_area=plot_area)
end

diametrictable(d::AbstractVector{<:Real}; plot_area::Union{Area,Real}=1.0) =
  diametrictable(d * DUNIT; plot_area=withunit(plot_area, AUNIT))

"""
    diametrictable(g::Symbol, d::Symbol, data::AbstractDataFrame; plot_area::Real=1.0)

Creates a diametric table for grouped data in a DataFrame.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `d::S`: The symbol representing the diameter values variable.
- `data::AbstractDataFrame`: The DataFrame containing the data.
- `plot_area::Area=1.0u"ha"`: The sampled plot area. A plain number is taken to be hectares.

# Units

The diameter column may hold either plain numbers (taken to be centimeters) or `Unitful`
quantities. The output columns carry units in both cases.

# Returns
- `DataFrame`: A DataFrame containing the diametric table for each group.

# Example
```julia-repl
julia> using DataFrames

# Define the data'
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]u"cm";
julia> species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"];
julia> data = DataFrame(species=species, diameters=diameters);

# Calculate the diametric table by group with with auto class width and plot area of 0.05 ha (500 m²)
julia> diametrictable(:species, :diameters, data, plot_area=0.05u"ha")
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
function diametrictable(g::Symbol, d::Symbol, data::AbstractDataFrame; plot_area::Union{Area,Real}=1.0)
  combine(groupby(data, g)) do df
    diametrictable(df[:, d], plot_area=plot_area)
  end
end

"""
    diametrictable(g::Symbol, d::Symbol, hi::Union{Len,Real}, data::AbstractDataFrame; plot_area::Union{Area,Real}=1.0u"ha")

Creates a diametric table for grouped data in a DataFrame with a specified class width.

# Arguments
- `g::S`: The symbol representing the grouping variable.
- `d::S`: The symbol representing the diameter values variable.
- `hi::Len`: The class width. A plain number is taken to be centimeters, or the unit of the diameter column when it carries one.
- `data::AbstractDataFrame`: The DataFrame containing the data.
- `plot_area::Area=1.0u"ha"`: The sampled plot area. A plain number is taken to be hectares.

# Returns
- `DataFrame`: A DataFrame containing the diametric table for each group.

# Example
```julia-repl
julia> using DataFrames

# Define the data'
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];
julia> species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"];
julia> data = DataFrame(species=species, diameters=diameters);

# Calculate the diametric table by group with with The class width of 3 and plot area of 0.05 ha (500 m²)
julia> diametrictable(:species, :diameters, 3, data, plot_area=0.05)
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
function diametrictable(g::Symbol, d::Symbol, hi::Union{Len,Real}, data::AbstractDataFrame; plot_area::Union{Area,Real}=1.0)
  combine(groupby(data, g)) do df
    diametrictable(df[:, d], hi, plot_area=plot_area)
  end
end