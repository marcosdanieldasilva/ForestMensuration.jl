"""
    basalarea(d::Len)

calculates the individual cross-sectional area (g) of a tree stem.

# arguments
- `d`: diameter at breast height as a unitful quantity e.g. 10u"cm" or 12u"inch". The diameter must be a positive value.

# returns
- `Quantity`: area in square feet ft2 if input is imperial or square meters m2 otherwise.

# mathematical basis
the calculation uses the standard geometric formula
```math
g = \\frac{\\pi d^2}{4}
```

# accuracy note

most stems are not perfectly circular
using a diameter tape slightly overestimates the true area because the circle
is the geometric figure with the smallest perimeter for a given area.

# examples

```julia
julia> basalarea(30u"cm")
0.07068583470577035 m^2

julia> basalarea(11.8u"inch")
0.7594363907740327 ft^2
```
"""
function basalarea(d::Len)
  ustrip(d) <= 0 && throw(DomainError("The diameter must be a positive value, observed values is $d"))

  g = quartπ * abs2(d)

  if unit(d) isa ImperialUnits
    return uconvert(u"ft^2", g)
  else
    return uconvert(u"m^2", g)
  end
end

"""
    basalarea(d::Real)

Calculates the basal area (g) of a tree given its diameter in centimeters.

# Description
This function computes the basal area of a tree, which is the cross-sectional area of the tree trunk at breast height (usually measured at 1.3 meters above ground). Basal area is a critical parameter in forest mensuration, used for estimating stand density, timber volume, and assessing competition among trees in a forest stand.

# Arguments
- `d::Real`: The diameter at breast height of the tree in **centimeters**. The diameter must be a positive value.

# Returns
- `Float64`: The basal area of the tree in **square meters**.

# Example
```julia
# Calculate the basal area for a tree with a diameter of 30 cm
julia> basalarea(30)
0.07068583470577035 m^2
```
"""
basalarea(d::Real) = basalarea(d * u"cm")

"""
    dm(d::AbstractVector{<:Len})

calculates the arithmetic mean diameter dbar of a forest stand or sample.

# arguments
- `d`: vector of diameters at breast height with length units (e.g., cm or inch).

# returns
- `Quantity`: the arithmetic mean diameter in the same units as input.

# mathematical basis
the simple average of all measured diameters
```math
\\bar{d} = \\frac{\\sum_{i=1}^{n} d_i}{n}
```

# technical description

the arithmetic mean is primarily used to characterize the frequency distribution of a stand
it is essential for biological research and serves as the base for hohenadl diameters
in even aged stands the arithmetic mean is generally smaller than the quadratic mean diameter dg
following the sequence dminus < dm < dg < dw < dz < dplus
unlike dg the arithmetic mean is highly sensitive to thinning interventions

# examples

```julia
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> dm(diameters)
17.25 cm
```
"""
dm(d::AbstractVector{<:Len}) = mean(d)

dm(d::AbstractVector{<:Real}) = dm(d * u"cm")

"""
    dg(d::AbstractVector{<:Len})

Calculate the quadratic mean diameter (dg) of a forest stand or sample. This metric represents the diameter of the tree with the mean basal area.

# arguments
- `d`: vector of diameters at breast height with length units (e.g., cm or inch).

# returns
- `Quantity`: quadratic mean diameter in the same unit as input.

# mathematical basis
calculated as the square root of the average of squared diameters
```math
d_g = \\sqrt{\\frac{\\sum_{i=1}^{n} d_i^2}{n}}
```

alternatively derived from arithmetic mean and variance

```math
d_g = \\sqrt{\\bar{d}^2 + s_d^2}
```

# technical description

dg is fundamental in forest mensuration and production models
it is preferred over arithmetic mean because it correlates stronger with stand volume
it gives additional weight to larger diameters which constitute more biomass
in even aged stands dg is always greater than or equal to the arithmetic mean.

# examples

```julia
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> dg(diameters)
17.779904386694547 cm
```
"""
dg(d::AbstractVector{<:Len}) = √(mean(abs2, d))

dg(d::AbstractVector{<:Real}) = dg(d * u"cm")

"""
    dw(d::AbstractVector{<:Len})

calculates weises diameter dw corresponding to the 60th percentile of the cumulative frequency distribution.

# arguments
- `d`: vector of diameters at breast height with length units (e.g., cm or inch).

# returns
- `Quantity`: diameter value at the 60th percentile in the same unit as input.

# mathematical basis
obtained from the inverse cumulative distribution function at 0.6
```math
d_w = F^{-1}(0.6)
```

# technical description

proposed by Weise 1880 this diameter approximates the tree with mean stand volume
in even aged stands it separates the 60 percent smallest trees from the 40 percent largest.

# examples

```julia
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> dw(diameters)
18.6 cm
```
"""
dw(d::AbstractVector{<:Len}) = quantile(d, 0.6)

dw(d::AbstractVector{<:Real}) = dw(d * u"cm")

"""
    dz(d::AbstractVector{<:Len})

calculates the basal area central diameter dz also known as the diameter of the median basal area.

# arguments
- `d`: vector of diameters at breast height with length units (e.g., cm or inch).

# returns
- `Quantity`: central diameter value in the same unit as input.

# mathematical basis
represents the diameter corresponding to the median individual basal area
mathematically simplified to the square root of the median of squared diameters
```math
d_z = \\sqrt{\\text{median}(d^2)}

```

# technical description

the central diameter is defined as the value corresponding to the median of the basal area distribution
it is a key reference diameter in relascopy bitterlichs method
it is generally larger than the arithmetic mean diameter and often close to the quadratic mean diameter.

# examples

```julia
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> dz(diameters)
17.266296649832007 cm

```
"""
dz(d::AbstractVector{<:Len}) = abs2.(d) |> median |> sqrt

dz(d::AbstractVector{<:Real}) = dz(d * u"cm")


# Auxiliary function that calculates the target number of dominant trees (k) for the plot.
# Returns NaN if the sample size is insufficient.
function dominantTreeCount(d::AbstractVector{<:Len}, area::Area)
  # detect unit system
  u = unit(d[1])
  if u isa ImperialUnits
    # imperial standard 40 trees per acre
    ntree = round(Int, ustrip(uconvert(u"ac", 40 * area)))
  else
    # metric standard 100 trees per hectare
    ntree = round(Int, ustrip(uconvert(u"ha", 100 * area)))
  end
  # validate against sample size
  # return nan if we need more trees than available
  if ntree <= 0 || ntree > length(d)
    return NaN64
  end
  return ntree
end

"""
    dd(d::AbstractVector{<:Len}, area::Area)

calculates the dominant diameter ddom based on Assmann's definition
average of the 100 thickest trees per hectare or 40 per acre.

# arguments
- `d`: vector of diameters at breast height with length units (e.g., cm or inch).
- `area`: total sampled area used to determine the number of dominant trees.

# returns
- `Quantity`: the arithmetic mean of the dominant trees in the same unit as input.

# technical description
the dominant diameter is a stable stand variable used to characterize site productive capacity
it uses the spatial density standard of 100 trees per hectare for metric systems
or 40 trees per acre for imperial systems
if the calculated number of trees exceeds the sample size the mean of the entire sample is returned.

# mathematical basis
the number of trees k to select is determined by
metric k = 100 * area_ha
imperial k = 40 * area_ac
then
```math
d_{dom} = \\frac{\\sum_{i=1}^{k} d_{i}}{k}

```

where di are the diameters sorted in descending order

# examples

```julia
# metric example 100 per ha
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> plotArea = 500u"m^2";
julia> dd(diameters, area)
21.0 cm
```
"""
function dd(d::AbstractVector{<:Len}, area::Area)
  # get target count of dominant trees
  ntree = dominantTreeCount(d, area)
  # return nan if sample is insufficient or invalid
  isnan(ntree) && return NaN64
  # calculate mean of the ntree thickest trees
  return dm(partialsort(d, 1:Int(ntree), rev=true))
end

dd(d::AbstractVector{<:Real}, area::Real) = dd(d * u"cm", area * u"ha")

"""
    dh(d::AbstractVector{<:Len})

calculates hohenadls diameters dminus and dplus
these metrics represent one standard deviation below and above the mean diameter.

# arguments
- `d`: vector of diameters at breast height with length units (e.g., cm or inch).

# returns
- `NamedTuple`: containing keys `lower` (d-) and `upper` (d+) in the same units as input.

# mathematical basis
calculated using arithmetic mean dbar and standard deviation s
```math
d_- = \\bar{d} - s
d_+ = \\bar{d} + s

```

# technical description

lower hohenadl diameter dminus is calculated as dbar minus s where dbar is the mean diameter and s is the standard deviation
approximately 16 percent of the trees have diameters below dminus representing one standard deviation below the mean
upper hohenadl diameter dplus is calculated as dbar plus s
approximately 84 percent of the trees have diameters below dplus representing one standard deviation above the mean
these metrics are useful for understanding the variability and distribution of diameters within the stand

# examples

```julia
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> hohenadl = dh(diameters)                        
(dl = 12.708524468853764 cm, du = 21.791475531146236 cm)
# access specific values
julia> hohenadl.dl
12.708524468853764 cm
julia> hohenadl.du
21.791475531146236 cm

```
"""
function dh(d::AbstractVector{<:Len})
  d̅ = dm(d)
  s = std(d, mean=d̅)
  d₋, d₊ = d̅ - s, d̅ + s
  return (; dl=d₋, du=d₊)
end

dh(d::AbstractVector{<:Real}) = dh(d * u"cm")

"""
    dmetrics(d::AbstractVector{<:Len}, area::Area=0.0u"ha")

Calculate a comprehensive set of dendrometric averages for a forest stand.
Aggregates the arithmetic mean, quadratic mean, Hohenadl, Weise, central basal area, and dominant diameters.

# Arguments
- `d`: vector of diameters at breast height with length units (e.g., cm or inch).
- `area`: sampled area used for dominant diameter calculation. Defaults to 0.0ha (which yields NaN for dominant diameter).

# Returns
- `DataFrame`: a single row dataframe containing:
  - `dl`: lower Hohenadl diameter (mean minus standard deviation), representing the lower variability threshold (approx. 16% of trees are smaller).
  - `dm`: arithmetic mean diameter, a basic measure of central tendency sensitive to thinning.
  - `dg`: quadratic mean diameter, the diameter of the tree with mean basal area (closely correlated with volume).
  - `dw`: Weise's diameter (60th percentile), approximating the tree with mean volume.
  - `dz`: central basal area diameter, the median of the basal area distribution (key for relascopy).
  - `dd`: dominant diameter, the mean of the thickest trees (100 per ha or 40 per ac), characterizing site capacity.
  - `du`: upper Hohenadl diameter (mean plus standard deviation), representing the upper variability threshold (approx. 84% of trees are smaller).
  - `dv`: coefficient of variation of the diameters, expressed as a percentage, representing the relative dispersion of the stand.

# Examples
```julia
julia> dmetrics(diameters)
1×8 DataFrame
 Row │ dl          dm          dg          dw          dz          dd       du          dv      
     │ Quantity…   Quantity…  Quantity…   Quantity…   Quantity…   Float64  Quantity…   Float64 
─────┼────────────────────────────────────────────────────────────────────────────────────────
   1 │ 12.7085 cm   17.25 cm  17.7799 cm    18.6 cm  17.2663 cm      NaN  21.7915 cm  26.3274
```

"""
function dmetrics(d::AbstractVector{<:Len}, area::Area=0.0u"ha")
  d̅ = dm(d)
  s = std(d, mean=d̅)
  d₋, d₊ = d̅ - s, d̅ + s
  return DataFrame(
    dl=d₋,
    dm=d̅,
    dg=dg(d),
    dw=dw(d),
    dz=dz(d),
    dd=dd(d, area),
    du=d₊,
    dv=s / d̅ * 100,
  )
end


dmetrics(d::AbstractVector{<:Real}, area::Real=0.0) = dmetrics(d * u"cm", area * u"ha")

"""
    hm(h::AbstractVector{<:Len})

Calculate the arithmetic mean height of a forest stand or sample.

# Arguments
- `h`: vector of tree heights with length units (e.g., m or ft).

# Returns
- `Quantity`: the arithmetic mean height in the same units as the input.

# Mathematical basis
The simple average of all measured heights:
```math
\\bar{h} = \\frac{\\sum_{i=1}^{n} h_i}{n}

```

# Technical description

The arithmetic mean height is a basic descriptive statistic of the stand's vertical structure. While useful for general characterization, it is highly sensitive to the presence of suppressed trees and thinning from below. In uneven-aged stands, it often underestimates the visual canopy height.

# Examples

```julia
julia> heights = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2] * u"m";

julia> hm(heights)
15.55 m

```

"""
hm(h::AbstractVector{<:Len}) = mean(h)

hm(d::AbstractVector{<:Real}) = hm(d * u"m")

"""
    hd(d::AbstractVector{<:Len}, h::AbstractVector{<:Len}, area::Area)

Calculate the dominant height (hd) of a forest stand based on Assmann's definition.

# Arguments
* `d`: vector of diameters at breast height with length units (e.g., cm or inch).
* `h`: vector of tree heights with length units (e.g., m or ft).
* `area`: total sampled area used to determine the number of dominant trees.

# Returns
* `Quantity`: the arithmetic mean height of the dominant trees in the same units as the input `h`. Returns `NaN` if the sample is insufficient.

# Mathematical basis
The number of dominant trees () is defined by the spatial density standard (100 per hectare or 40 per acre). The dominant height is the mean height of those  thickest trees:

```math
h_{dom} = \\frac{\\sum_{i=1}^{k} h_{i}}{k}

```

where  corresponds to the heights of the trees with the largest diameters.

# Technical description

Dominant height is one of the most important variables in forest mensuration. Because it is practically unaffected by thinning interventions and stand density, it is universally used to evaluate site quality (Site Index) and to construct growth and yield models.

# Examples

```julia
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> heights = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2] * u"m";
julia> plotArea = 500u"m^2";
julia> hd(diameters, heights, plotArea)
18.5 m
```

"""
function hd(d::AbstractVector{<:Len}, h::AbstractVector{<:Len}, area::Area)
  # validate vector lengths
  length(h) == length(d) || throw(DimensionMismatch("vectors must have same length"))
  # get target count of dominant trees
  ntree = dominantTreeCount(d, area)
  # return nan if sample is insufficient or invalid
  isnan(ntree) && return NaN64
  # get indices of the ntree thickest trees
  idx = partialsortperm(d, 1:Int(ntree), rev=true)
  return hm(h[idx])
end

"""
    hg(d::AbstractVector{<:Len}, h::AbstractVector{<:Len})

Calculate Lorey's mean height (hg), which is the mean height weighted by the basal area of each tree.

# Arguments
* `d`: vector of diameters at breast height with length units (e.g., cm or inch).
* `h`: vector of tree heights with length units (e.g., m or ft).

# Returns
* `Quantity`: Lorey's mean height in the same units as the input `h`.

# Mathematical basis
Calculated by dividing the sum of the products of individual heights and basal areas by the total stand basal area:

```math
h_g = \\frac{\\sum_{i=1}^{n} g_i h_i}{\\sum_{i=1}^{n} g_i}

```

# Technical description
Proposed by Lorey (1878), this is the standard mean height used in volumetric calculations and forest production models. By weighting heights according to basal area, it reduces the influence of smaller, suppressed trees and emphasizes the larger trees that constitute the bulk of the stand's volume. It is invariably larger than the arithmetic mean height.

# Examples

```julia
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> heights = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2] * u"m";
julia> hg(heights, diameters)
17.148042704626334 m
```

"""
function hg(d::AbstractVector{<:Len}, h::AbstractVector{<:Len})
  if length(h) != length(d)
    throw(DimensionMismatch("height and diameter vectors must have the same length"))
  end
  g = basalarea.(d)
  return sum(h .* g) / sum(g)
end

"""
    hmetrics(d::AbstractVector{<:Len}, h::AbstractVector{<:Len}, area::Area=0.0u"ha")

Calculate a comprehensive set of vertical structure metrics for a forest stand.
Aggregates the lower boundary, arithmetic mean, dominant height, Lorey's mean height, upper boundary, and relative dispersion.

# Arguments
* `d`: vector of diameters at breast height with length units (e.g., cm or inch).
* `h`: vector of tree heights with length units (e.g., m or ft).
- `area`: sampled area used for dominant diameter calculation. Defaults to 0.0ha (which yields NaN for dominant diameter).

# Returns
* `DataFrame`: a single row dataframe containing:
* `hl`: lower height boundary (mean minus standard deviation).
* `hm`: arithmetic mean height.
* `hd`: dominant height (Assmann's definition, based on thickest trees).
* `hg`: Lorey's mean height (weighted by basal area).
* `hu`: upper height boundary (mean plus standard deviation).
* `hv`: coefficient of variation of the heights, expressed as a percentage.

```julia
julia> heights = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2] * u"m";
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> plotArea = 500u"m^2";
julia> hmetrics(heights, diameters, plotArea)
1×6 DataFrame
 Row │ hl         hm         hd         hg         hu         hv      
     │ Quantity…  Quantity…  Quantity…  Quantity…  Quantity…  Float64
─────┼────────────────────────────────────────────────────────────────
   1 │ 11.9589 m    15.55 m     18.5 m   17.148 m  19.1411 m   23.094
```
"""
function hmetrics(d::AbstractVector{<:Len}, h::AbstractVector{<:Len}, area::Area=0.0u"ha")
  # validate vector lengths
  length(h) == length(d) || throw(DimensionMismatch("height and diameter vectors must have same length"))
  # mean height h̅
  h̅ = hm(h)
  s = std(h, mean=h̅)
  h₋, h₊ = h̅ - s, h̅ + s
  return DataFrame(
    hl=h₋, hm=h̅,
    hd=hd(h, d, area),
    hg=hg(h, d),
    hu=h₊,
    hv=s / h̅ * 100,
  )
end

"""
standmetrics(d::AbstractVector{<:Len}, h::AbstractVector{<:Len}, area::Area)

Calculate the complete summary of inventory metrics for a forest plot, scaling values to per-hectare or per-acre equivalents.

# Arguments
* `d`: vector of diameters at breast height with length units (e.g., cm or inch).
* `h`: vector of tree heights with length units (e.g., m or ft).
- `area`: sampled area used for dominant diameter calculation. Defaults to 0.0ha (which yields NaN for dominant diameter).

# Returns
* `DataFrame`: a single row dataframe aggregating absolute plot totals, scaled per-area values, and all diameter and height metrics. Contains:
* `n`: total number of sampled trees in the plot.
* `g`: total basal area of the sampled plot.
* `EF`: expansion factor used to scale plot data to a per-hectare or per-acre basis.
* `N`: extrapolated number of trees per hectare/acre.
* `G`: extrapolated total basal area per hectare/acre.
* Includes all metrics from `dmetrics` (dm, dg, dd, etc.) and `hmetrics` (hm, hg, hd, etc.).


# Technical description

This is the primary entry point for characterizing an individual plot in a forest inventory. It unifies horizontal structure (diameters), vertical structure (heights), and spatial density (extrapolated totals). The Expansion Factor (EF) automatically detects whether the inputs are in the metric or imperial system and scales accordingly.

# Examples

```julia
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm";
julia> heights = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2] * u"m";
julia> plotArea = 500u"m^2";
julia> standmetrics(diameters, heights, plotArea)
1×19 DataFrame
 Row │ n      g             EF          N            G                  dl          dm         dg          dw         dz          dd         du          dv       hl          hm         hd         hg         hu          hv      
     │ Int64  Quantity…     Quantity…   Quantity…    Quantity…          Quantity…   Quantity…  Quantity…   Quantity…  Quantity…   Quantity…  Quantity…   Float64  Quantity…   Quantity…  Quantity…  Quantity…  Quantity…   Float64 
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │    10  0.248284 m^2  20.0 ha^-1  200.0 ha^-1  4.96568 m^2 ha^-1  12.7085 cm   17.25 cm  17.7799 cm    18.6 cm  17.2663 cm    21.0 cm  21.7915 cm  26.3274  12.7085 cm   17.25 cm     18.5 m   17.148 m  21.7915 cm  26.3274
```
"""
function standmetrics(d::AbstractVector{<:Len}, h::AbstractVector{<:Len}, area::Area)
  n = length(d)
  g = basalarea.(d) |> sum
  EF = unit(g) == u"m^2" ? uconvert(u"ha", area)^-1.0 : uconvert(u"ac", area)^-1.0
  N = n * EF
  G = g * EF
  dmet = dmetrics(d, area)
  hmet = hmetrics(h, d, area)
  return hcat(
    DataFrame(
      n=n,
      g=g,
      EF=EF,
      N=N,
      G=G
    ),
    dmet,
    hmet
  )
end

standmetrics(d::AbstractVector{<:Real}, h::AbstractVector{<:Real}, area::Real) = standmetrics(d * u"cm", h * u"m", area * u"ha")
