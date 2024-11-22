# Partial function to calculate deltas for age and index age
# Calculate the delta (difference) between the model matrices for age and index age.
# Arguments
# - `model::StatsModels.TableRegressionModel`: The fitted regression model.
# - `data_age::AbstractDataFrame`: The data frame containing the current age data.
# - `index_age::Real`: The index age for which the delta is calculated.
# Returns
# - `Δ::Vector{Real}`: The delta for each observation.
function _calculate_delta(model::TableRegressionModel, data_age::AbstractDataFrame, index_age::Real)
  (yname, xname, qname...) = propertynames(model.mf.data)
  data_index_age = deepcopy(data_age[!, [yname, xname, qname...]])
  data_index_age[!, xname] .= index_age

  mm_age = modelmatrix(model.mf.f.rhs.terms[2:end], data_age)
  mm_index_age = modelmatrix(model.mf.f.rhs.terms[2:end], data_index_age)

  β = coef(model)[2:end]'

  # Calculate the difference in model matrices
  Δ = mm_index_age .- mm_age

  # Compute the sum for each observation
  Δ = sum(β .* Δ, dims=2)[:]

  return Δ
end

"""
    site_classification(model::TableRegressionModel, data_age::AbstractDataFrame, index_age::Real)
    
Calculate the site classification (site index) for each observation given a fitted model, a table do predict and an index age.

# Description:

Site classification is a method used in forestry to evaluate the productivity and quality of a forest site. It is based on the relationship between the dominant height of trees and their age. By comparing the current dominant height and age of trees to the expected height at a standard reference age (known as the index age), we can classify the site into different productivity classes.

- **Dominant Height**: The average height of the tallest trees in a stand, considered indicators of site productivity.
- **Age**: The current age of the trees in the stand.
- **Index Age**: A standard age used for comparison, often corresponding to the rotation age.

By predicting the dominant height at the index age using a fitted growth model, we obtain the site index. This index reflects the potential height growth of trees on that site and is a key parameter for forest management, aiding in decision-making regarding thinning, harvesting, and sustainability practices.

# Arguments:

- `model::StatsModels.TableRegressionModel`: The fitted regression model that relates dominant height to age. This model is typically derived from empirical data and captures the growth patterns of the species in question.
- `data_age::AbstractDataFrame`: A data frame containing the current age data for each observation. It should include an `age` column representing the age of the trees.
- `index_age::Real`: The index age (reference or rotation age) used for site classification. This is the age at which the site index is evaluated and compared.

# Returns
- `site::Vector{Real}`: A vector containing the site classification (site index) for each observation. Each element represents the expected dominant height at the index age, providing a measure of site productivity.


# Exemple
```
julia> using DataFrames
using DataFrames

# Create a DataFrame containing tree plot data
julia> data = DataFrame(
    plot = repeat(1:6, inner=5),
    age = repeat([36, 48, 60, 72, 84], outer=6),
    h = [13.6, 17.8, 21.5, 21.5, 21.8,
          14.3, 17.8, 21.0, 21.0, 21.4,
          14.0, 17.5, 21.2, 21.2, 21.4,
          13.4, 18.0, 20.8, 20.8, 23.2,
          13.2, 17.4, 20.3, 20.3, 22.0,
          13.2, 17.8, 21.3, 21.3, 22.5]
);

# Perform regression to model the relationship between height (`h`) and age (`age`)
# `criteria_selection` selects the best regression model based on predefined criteria
julia> reg = regression(:h, :age, data) |> criteria_selection;

# Define the data set to predict the site classification
julia> data_to_predict = DataFrame(
    plot = repeat(7:8, inner=5),
    age = repeat([36, 48, 60, 72, 84], outer=2),
    h = [13.1, 17.4, 20.8, 20.8, 22.6,
          14.3, 17.3, 20.4, 20.4, 20.9]
);

# Use the selected regression model to classify site quality for each plot
# at an age of 60 months (5 years). Returns a vector of predicted site classification.
julia> site_classification(reg, data_to_predict, 60)
10-element Vector{Float64}:
 19.3
 19.7
 20.8
 19.8
 21.1
 22.2
 19.6
 20.4
 19.4
 19.6
)
```
"""
function site_classification(model::TableRegressionModel, data_age::AbstractDataFrame, index_age::Real)
  if index_age <= 0
    throw(DomainError("Index Age must be positive."))
  end
  Δ = _calculate_delta(model, data_age, index_age)
  # Calculate the site classification
  y = model.mf.f.lhs
  hdom = modelcols(y, data_age)
  site = @. hdom + Δ

  if isa(y, FunctionTerm)
    site = _prediction(model, [index_age], site)
  end

  return round.(site, digits=1)
end

"""
    site_classification(model::TableRegressionModel, index_age::Real)

Calculate the site classification (site index) for each observation given a fitted model and an index age.

# Description:

Site classification is a method used in forestry to evaluate the productivity and quality of a forest site. It is based on the relationship between the dominant height of trees and their age. By comparing the current dominant height and age of trees to the expected height at a standard reference age (known as the index age), we can classify the site into different productivity classes.

- **Dominant Height**: The average height of the tallest trees in a stand, considered indicators of site productivity.
- **Age**: The current age of the trees in the stand.
- **Index Age**: A standard age used for comparison, often corresponding to the rotation age.

By predicting the dominant height at the index age using a fitted growth model, we obtain the site index. This index reflects the potential height growth of trees on that site and is a key parameter for forest management, aiding in decision-making regarding thinning, harvesting, and sustainability practices.

# Arguments:

- `model::StatsModels.TableRegressionModel`: The fitted regression model that relates dominant height to age. This model is typically derived from empirical data and captures the growth patterns of the species in question.
- `index_age::Real`: The index age (reference or rotation age) used for site classification. This is the age at which the site index is evaluated and compared.

# Returns
- `site::Vector{Real}`: A vector containing the site classification (site index) for each observation. Each element represents the expected dominant height at the index age, providing a measure of site productivity.

# Exemple
```
julia> using DataFrames
using DataFrames

# Create a DataFrame containing tree plot data
julia> data = DataFrame(
    plot = repeat(1:6, inner=5),
    age = repeat([36, 48, 60, 72, 84], outer=6),
    h = [13.6, 17.8, 21.5, 21.5, 21.8,
          14.3, 17.8, 21.0, 21.0, 21.4,
          14.0, 17.5, 21.2, 21.2, 21.4,
          13.4, 18.0, 20.8, 20.8, 23.2,
          13.2, 17.4, 20.3, 20.3, 22.0,
          13.2, 17.8, 21.3, 21.3, 22.5]
);

# Perform regression to model the relationship between height (`h`) and age (`age`)
# `criteria_selection` selects the best regression model based on predefined criteria
julia> reg = regression(:h, :age, data) |> criteria_selection;

# Use the selected regression model to classify site quality for each plot
# at an age of 60 months (5 years). Returns a vector of predicted site classification.
julia> site_classification(reg, 60)
30-element Vector{Float64}:
 20.5
 20.3
 21.5
 20.4
 20.4
 22.2
 20.3
 21.0
 19.9
 20.0
 21.4
 19.9
 21.2
 20.1
 20.0
 20.0
 20.5
 20.8
 19.8
 21.6
 19.5
 19.7
 20.3
 19.3
 20.5
 19.5
 20.3
 21.3
 20.2
 21.0
)
```
"""
function site_classification(model::TableRegressionModel, index_age::Real)
  # Extract the data from the fitted model
  data_age = model.mf.data |> DataFrame
  return site_classification(model, data_age, index_age)
end

"""
    hdom_classification(model::TableRegressionModel, data_age::AbstractDataFrame, index_age::Real, site::Vector{<:Real})

Calculate the dominant height for each observation given the site classification, a fitted model, and an index age.

# Description:

The `hdom_classification` function calculates the dominant height for each observation based on the provided site classification values, using a fitted growth model and the index age. This function essentially works in the reverse direction of site classification; given the site index (site classification) and the ages, it predicts the dominant heights. It is useful for forecasting tree growth and estimating future stand development based on site productivity classes.

- **Dominant Height**: The expected average height of the tallest trees in a stand at a given age.
- **Age**: The current age of the trees for which the dominant height is to be calculated.
- **Index Age**: A standard age used as a reference point in the growth model.
- **Site Classification (Site Index)**: Values representing the site productivity class, typically the expected dominant height at the index age.

# Arguments:

- `model::StatsModels.TableRegressionModel`: The fitted regression model relating dominant height to age and site index. This model is used to predict heights based on site classes.
- `data_age::AbstractDataFrame`: A DataFrame containing the current age data. It should include an `age` column representing the ages for which the dominant height is to be calculated.
- `index_age::Real`: The index age used in the growth model. It serves as the reference age for site index calculations.
- `site::Vector{<:Real}`: A vector containing the site classification values (site indices) for each observation. Each value corresponds to an expected dominant height at the index age.

# Returns:

- `hdom::Vector{Real}`: A vector containing the predicted dominant heights for each observation at their respective ages.

# Exemple
```
julia> using DataFrames
using DataFrames

# Create a DataFrame containing tree plot data
julia> data = DataFrame(
    plot = repeat(1:6, inner=5),
    age = repeat([36, 48, 60, 72, 84], outer=6),
    h = [13.6, 17.8, 21.5, 21.5, 21.8,
          14.3, 17.8, 21.0, 21.0, 21.4,
          14.0, 17.5, 21.2, 21.2, 21.4,
          13.4, 18.0, 20.8, 20.8, 23.2,
          13.2, 17.4, 20.3, 20.3, 22.0,
          13.2, 17.8, 21.3, 21.3, 22.5]
);

# Perform regression to model the relationship between height (`h`) and age (`age`)
# `criteria_selection` selects the best regression model based on predefined criteria
julia> reg = regression(:h, :age, data) |> criteria_selection;

# Define the data set to predict the site classification
julia> data_to_predict = DataFrame(
    plot = repeat(7:8, inner=5),
    age = repeat([36, 48, 60, 72, 84], outer=2),
    h = [13.1, 17.4, 20.8, 20.8, 22.6,
          14.3, 17.3, 20.4, 20.4, 20.9]
);

# Use the selected regression model to classify site quality for each plot
# at an age of 60 months (5 years). Returns a vector of predicted site classification.
julia> sites = site_classification(reg, data_to_predict, 60);

# The predicted dominant height for each observation
julia> hdom_classification(reg, data_to_predict, 60, sites)
10-element Vector{Float64}:
 13.1
 17.4
 20.8
 20.9
 22.7
 14.3
 17.3
 20.4
 20.4
 20.9
```
"""
function hdom_classification(model::TableRegressionModel, data_age::AbstractDataFrame, index_age::Real, site::Vector{<:Real})
  if index_age <= 0
    throw(DomainError("Index Age must be positive."))
  elseif any(x -> x < 0, site)
    throw(DomainError("Site values must be positive"))
  end
  Δ = _calculate_delta(model, data_age, index_age)
  (yname, xname, qname...) = propertynames(model.mf.data)
  site_data = deepcopy(data_age[!, [yname, xname, qname...]])
  site_data[!, yname] .= site
  site_data[!, xname] .= index_age
  # Calculate the site classification
  y = model.mf.f.lhs
  site = modelcols(y, site_data)
  hdom = @. site - Δ

  if isa(y, FunctionTerm)
    hdom = _prediction(model, data_age[!, xname], hdom)
  end

  return round.(hdom, digits=1)
end

"""
    site_table(model::TableRegressionModel, index_age::Real)
    site_table(model::TableRegressionModel, index_age::Real, hi::Real)

Calculate the site table and plot a site graph given a fitted model, index age, and height increment.

# Description:

The `site_table` function generates a site table based on a fitted growth model, an index age, and a specified height increment. A site table is a tool used in forestry to estimate the expected dominant height of trees at various ages for different site quality classes. It provides a tabulated representation of height growth over time, allowing foresters to assess site productivity and make informed management decisions.

- **Site Quality Classes**: These are categories that represent different levels of site productivity, typically based on the site index (the expected dominant height at the index age).
- **Index Age**: A standard reference age used for site classification and comparison between different sites.
- **Height Increment (`hi`)**: The increment by which site index classes are divided. It determines the range and granularity of the site quality classes in the table.

By using the fitted model to predict dominant heights at various ages and site indices, the site table provides a comprehensive overview of growth patterns across different site qualities.

# Arguments
- `model::StatsModels.TableRegressionModel`: The fitted regression model.
- `index_age::Real`: The index age for site classification.
- `hi::Real`: The height increment for site classification.

# Returns
- `SiteAnalysis`: A struct containing a dataFrame containing the site table with predicted dominant heights at various ages for different site index classes and a site plot.

# Exemple
```
julia> using DataFrames
using DataFrames

# Create a DataFrame containing tree plot data
julia> data = DataFrame(
    plot = repeat(1:6, inner=5),
    age = repeat([36, 48, 60, 72, 84], outer=6),
    h = [13.6, 17.8, 21.5, 21.5, 21.8,
          14.3, 17.8, 21.0, 21.0, 21.4,
          14.0, 17.5, 21.2, 21.2, 21.4,
          13.4, 18.0, 20.8, 20.8, 23.2,
          13.2, 17.4, 20.3, 20.3, 22.0,
          13.2, 17.8, 21.3, 21.3, 22.5]
);

# Perform regression to model the relationship between height (`h`) and age (`age`)
# `criteria_selection` selects the best regression model based on predefined criteria
julia> reg = regression(:h, :age, data) |> criteria_selection;

# Use the selected regression model to generate the site table and site plot
julia> site_table(reg, 60, 1)
5×5 DataFrame
 Row │ age      S_19.5   S_20.5   S_21.5   S_22.5  
     │ Float64  Float64  Float64  Float64  Float64
─────┼─────────────────────────────────────────────
   1 │    36.0     13.2     13.6     14.0     14.4
   2 │    48.0     17.2     18.0     18.7     19.5
   3 │    60.0     19.5     20.5     21.5     22.5
   4 │    72.0     20.5     21.6     22.8     23.9
   5 │    84.0     20.8     22.0     23.1     24.3
```
"""
function site_table(model::TableRegressionModel, index_age::Real, hi::Real)
  if hi <= 0
    throw(DomainError("The height increment must be positive."))
  end
  # Extract property names from the fitted model's data
  (hd, age, q...) = propertynames(model.mf.data)
  # Calculate site classification
  site = site_classification(model, index_age)
  # Get unique and sorted site classes
  sites = _class_center.(site, hi) |> unique |> sort
  # Get unique and sorted ages
  ages = model.mf.data[age] |> unique |> sort
  # Repeat ages for each site class
  repeated_ages = repeat(ages, outer=length(sites))
  # Create repeated site classes
  repeated_sites = vcat([fill(s, length(ages)) for s in sites]...)
  # Create initial DataFrame with repeated ages and sites
  site_table = DataFrame(age => repeated_ages, hd => repeated_sites)
  # Predict the dominant heights for the site classification
  hdom_predict = hdom_classification(model, site_table, index_age, repeated_sites)
  # Insert the predicted heights into the DataFrame
  insertcols!(site_table, :site => hdom_predict, makeunique=true)
  # transform into a pivote table by site
  site_table = unstack(site_table, propertynames(site_table)...) |> dropmissing!
  # Rename columns to reflect site classes
  new_column_names = [Symbol("S_$(s)") for s in names(site_table)[2:end]]
  rename!(site_table, [age; new_column_names])
  # Create the site plot
  site_plot = Plots.plot(
    repeated_ages, hdom_predict,
    group=categorical(repeated_sites, levels=sort(sites, rev=true)),
    xlabel=age,
    ylabel=hd,
    title="Site Classification Plot",
    legend=:outertopright,
    markerstrokewidth=0,
    seriesalpha=0.6,
    framestyle=:box,
    margin=3mm,
    color_palette=cgrad(:darktest, categorical=true),
    tick_direction=:out,
    grid=:none,
    titlefont=10,
    fontfamily="times",
    guidefontsize=9,
    legendfontsize=7
  )

  return SiteAnalysis(site_table, site_plot)

end

function site_table(model::TableRegressionModel, index_age::Real)
  # Calculate site classification
  site = site_classification(model, index_age)
  # Compute the class breadth
  h = _amplitude(site)
  k = _sturges(length(site))
  hi = _class_breadth(h, k)
  return site_table(model, index_age, hi)
end