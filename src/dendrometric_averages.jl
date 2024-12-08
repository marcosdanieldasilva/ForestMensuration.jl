"""
    basal_area(d::Real)

Calculates the basal area (g) of a tree given its diameter in centimeters.

# Description
This function computes the basal area of a tree, which is the cross-sectional area of the tree trunk at breast height (usually measured at 1.3 meters above ground). Basal area is a critical parameter in forest mensuration, used for estimating stand density, timber volume, and assessing competition among trees in a forest stand.

# Arguments
- `d::Real`: The diameter at breast height (DBH) of the tree in **centimeters**. The diameter must be a positive value.

# Returns
- `Float64`: The basal area of the tree in **square meters**.

# Example
```julia
# Calculate the basal area for a tree with a diameter of 30 cm
julia> basal_area(30.0)
0.07068583470577035
```
"""
@inline function basal_area(d::Real)
  if d <= 0
    throw(DomainError("Diameter must be positive."))
  else
    return (π / 40000) * d^2
  end
end

"""
    dendrometric_averages(d::Vector; area::Real=1.0)

Calculates various dendrometric averages of a forest stand, including mean diameter, quadratic mean 
diameter, Hohenadl's diameters, Weise's diameter, diameter of the tree with central basal area, and 
the mean diameter of the 100 largest trees per hectare.

# Description
This function computes several dendrometric averages based on a vector of tree diameters. These metrics 
  are essential for forest inventory and management, providing insights into the stand structure, volume 
    estimation, and growth patterns.

# Arguments
- `d::Vector{<:Real}`: A vector of diameters at breast height (DBH) of the trees in centimeters. All diameters must be positive values.
- `area::Real=1.0`: The area in hectares over which the diameters were measured. Default is 1.0 hectare.

# Returns
- `DataFrame`: A DataFrame containing the calculated dendrometric averages:
  - `d₋`: Lower Hohenadl's diameter, calculated as (d̄ - s), where d̄ is the mean diameter and s is the standard deviation. Approximately 16% of the trees have diameters below `d₋`. This metric represents one standard deviation below the mean diameter and is useful for understanding the variability and distribution of diameters within the stand.
  - `d̄`: Mean diameter, the arithmetic mean of the diameters. It is a basic measure of central tendency but can be influenced by thinning practices. The mean diameter is a fundamental measure but can be significantly affected by the removal of smaller or larger trees during thinning operations.
  - `dg`: Quadratic mean diameter, calculated using the mean basal area as `dg = sqrt((40000 * mean(g)) / π)`, where `g` is the basal area of each tree. It closely approximates the diameter of the tree with mean basal area and is less affected by extreme values or thinning. This metric provides a better estimate for volume calculations and is less sensitive to variations in the data.
  - `dw`: Weise's diameter, the diameter at which 60% of the trees have smaller diameters (the 60th percentile). It approximates the diameter of the tree with mean volume. Weise's diameter is considered a good approximation of the diameter of the tree with mean volume and is less influenced by thinning practices.
  - `dz`: Diameter of the tree with central basal area, calculated from the median basal area as `dz = sqrt((40000 * median(g)) / π)`. It represents the diameter that divides the stand's total basal area into two equal parts, effectively splitting the cumulative basal area. This metric is less influenced by the removal of smaller trees and provides insight into the stand's structure.
  - `d₁₀₀`: Mean diameter of the 100 largest trees per hectare. If there are fewer than 100 trees per hectare, it returns `NaN`. This metric provides insight into the size of the largest trees in the stand, which can be important for management objectives like timber production.
  - `d₊`: Upper Hohenadl's diameter, calculated as (d̄ + s). Approximately 84% of the trees have diameters below `d₊`. This metric represents one standard deviation above the mean diameter and is useful for understanding the variability and distribution of diameters within the stand.

# Example
```julia
julia> using DataFrames

# Sample diameters in centimeters
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];

# Calculate dendrometric averages
julia> dendrometric_averages(diameters, area=0.05)
1×7 DataFrame
 Row │ d₋       d̅        dg       dw       dz       d₁₀₀     d₊
     │ Float64  Float64  Float64  Float64  Float64  Float64  Float64
─────┼───────────────────────────────────────────────────────────────
   1 │ 12.7085    17.25  17.7799     18.6  17.2663     21.0  21.7915
```
"""
function dendrometric_averages(d::Vector; area::Real=1.0)
  if any(x -> x <= 0, d)
    throw(DomainError("All elements in vector d must be positive."))
  elseif area <= 0.0
    throw(DomainError("Area must be positive."))
  else
    d̅, s = mean_and_std(d)
    d₋, d₊ = d̅ - s, d̅ + s
    g = basal_area.(d)
    dg = √((40000 * mean(g)) / π)
    dw = quantile(d, 0.6)
    dz = √((40000 * median(g)) / π)
    n_tree = round(Int, (100 * area))
    d₁₀₀ = n_tree < length(d) ? mean(partialsort(d, 1:n_tree, rev=true)) : NaN64
    return DataFrame(d₋=d₋, d̅=d̅, dg=dg, dw=dw, dz=dz, d₁₀₀=d₁₀₀, d₊=d₊)
  end
end

"""
    dendrometric_averages(p::Symbol, d::Symbol, data::AbstractDataFrame; area::Real=1.0)

Calculates various dendrometric averages for each group in a dataset, grouping by a specified column, 
  and using the diameters specified in another column.

# Description
This function computes several dendrometric averages based on a specified diameter column in a DataFrame, 
  grouped by another specified column. These metrics are essential for forest inventory and management, 
  providing insights into the stand structure, volume estimation, and growth patterns for each group.

# Arguments
- `p::Symbol`: The symbol representing the column name in `data` used to group the data.
- `d::Symbol`: The symbol representing the column name in `data` that contains the diameters at breast height (DBH) of the trees in centimeters. All diameters must be positive values.
- `data::AbstractDataFrame`: The DataFrame containing the dataset with at least the columns specified by `p` and `d`.
- `area::Real=1.0`: The area in hectares over which the diameters were measured. Default is 1.0 hectare.

# Returns
- `DataFrame`: A DataFrame containing the calculated dendrometric averages for each group defined by `p`.

# Example
```julia
julia> using DataFrames

# Sample data
julia> data = DataFrame(
  species = ["Oak", "Oak", "Oak", "Pine", "Pine", "Pine"],
  diameter = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0]
);

# Calculate dendrometric averages grouped by species
julia> dendrometric_averages(:species, :diameter, data; area=0.05)
2×8 DataFrame
 Row │ species  d₋       d̅        dg       dw       dz       d₁₀₀     d₊
     │ String   Float64  Float64  Float64  Float64  Float64  Float64  Float64
─────┼────────────────────────────────────────────────────────────────────────
   1 │ Oak         10.5     12.0  12.0623     12.3     12.0      NaN     13.5
   2 │ Pine        15.0     16.5  16.5454     16.8     16.5      NaN     18.0
```
"""
function dendrometric_averages(p::Symbol, d::Symbol, data::AbstractDataFrame; area::Real=1.0)
  return combine(groupby(data, p)) do df
    dendrometric_averages(df[:, d], area=area)
  end
end

"""
    dendrometric_averages(model::FittedLinearModel; area::Real=1.0)
  
Calculates various dendrometric averages of a forest stand and estimates the corresponding heights for 
  each diameter using a regression model.

# Description
This function computes several dendrometric averages based on a regression model and estimates the heights
   associated with each calculated diameter using the provided regression model. These metrics are essential
    for forest inventory and management, providing insights into the stand structure, volume estimation, 
      and growth patterns.

# Arguments
- `model::FittedLinearModel`: A regression model used to predict heights from diameters. The model should be trained with diameters as predictors and heights as the response variable.
- `area::Real=1.0`: The area in hectares over which the diameters were measured. Default is 1.0 hectare.

# Returns
- `DataFrame`: A DataFrame containing the calculated dendrometric averages and the estimated heights. Heights are estimated by applying the regression model to each calculated diameter.
  - `h₋`: Estimated height corresponding to `d₋`.
  - `h̄`: Estimated height corresponding to `d̄`.
  - `hg`: Estimated height corresponding to `dg`.
  - `hw`: Estimated height corresponding to `dw`,.
  - `hz`: Estimated height corresponding to `dz`.
  - `d₁₀₀`: Estimated height corresponding to `d₁₀₀`.
  - `h₊`: Estimated height corresponding to `d₊`.

# Example
```julia
julia> using DataFrames

# Assume we have a trained regression model `height_model` that predicts height from diameter
# Sample data for the regression model
julia> heights = [8.0, 9.5, 11.0, 12.5, 14.0, 15.5, 17.0, 18.5, 20.0, 21.5];
julia> diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0];
julia> data = DataFrame(h = heights, d = diameters);

# Fit a linear regression model
julia> height_models = regression(:h, :d, data);
julia> best_model = criteria_selection(height_models);

# Calculate dendrometric averages and estimate heights
julia> dendrometric_averages(best_model; area=0.05)
1×14 DataFrame
 Row │ d₋       d̅        dg       dw       dz       d₁₀₀     d₊       h₋       h̅        hg       hw       hz       h₊       h₁₀₀
     │ Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 12.7085    17.25  17.7799     18.6  17.2663     21.0  21.7915  10.2085    14.75  15.2799     16.1  14.7663     18.5  19.2915
```
"""
function dendrometric_averages(model::FittedLinearModel; area::Real=1.0)
  data = model.data
  var_names = propertynames(data) |> collect
  if length(data) == 2
    result_table = dendrometric_averages(data[2], area=area)
    average_heights = map(i -> predict(model, DataFrame(var_names[2] => result_table[:, i]))[1], 1:7)
    return hcat(
      result_table,
      DataFrame(
        h₋=average_heights[1], h̅=average_heights[2], hg=average_heights[3],
        hw=average_heights[4], hz=average_heights[5], h₊=average_heights[6],
        h₁₀₀=average_heights[7]
      )
    )
  else
    result_table = combine(groupby(DataFrame(data), var_names[3:end])) do df
      dendrometric_averages(df[:, var_names[2]], area=area)
    end
    s = size(result_table, 2)
    average_heights = hcat(
      map(i -> predict(model, hcat(result_table[:, 1:s-7], DataFrame(var_names[2] => result_table[:, s-i]))), 6:-1:0)...
    )
    return hcat(result_table, DataFrame(average_heights, [:h₋, :h̅, :hg, :hw, :hz, :h₁₀₀, :h₊]))
  end
end