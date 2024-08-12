# Partial function to calculate deltas for age and index age
"""
Calculate the delta (difference) between the model matrices for age and index age.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted regression model.
- `data_age::AbstractDataFrame`: The data frame containing the current age data.
- `index_age::Real`: The index age for which the delta is calculated.

# Returns
- `Δ::Vector{Real}`: The delta for each observation.
"""
function _calculate_delta(fitted_model::FittedLinearModel, data_age::AbstractDataFrame, index_age::Real)
  (yname, xname, qname...) = propertynames(fitted_model.data)
  data_index_age = deepcopy(data_age[!, [yname, xname, qname...]])
  data_index_age[!, xname] .= index_age

  mm_age = modelmatrix(fitted_model.formula.rhs.terms[2:end], data_age)
  mm_index_age = modelmatrix(fitted_model.formula.rhs.terms[2:end], data_index_age)

  β = coef(fitted_model)[2:end]'

  # Calculate the difference in model matrices
  Δ = mm_index_age .- mm_age

  # Compute the sum for each observation
  Δ = sum(β .* Δ, dims=2)[:]
  return Δ
end

"""
Calculate the site classification given the fitted model and index age.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted regression model.
- `data_age::AbstractDataFrame`: The data frame containing the current age data.
- `index_age::Real`: The index age for site classification.

# Returns
- `site::Vector{Real}`: The site classification for each observation.
"""
function site_classification(fitted_model::FittedLinearModel, data_age::AbstractDataFrame, index_age::Real)
  Δ = _calculate_delta(fitted_model, data_age, index_age)
  # Calculate the site classification
  y = fitted_model.formula.lhs
  hdom = modelcols(y, data_age)
  site = hdom .+ Δ
  if isa(y, FunctionTerm)
    site = _predict(site, [index_age], fitted_model.σ², nameof(y.f))
  end

  return round.(site, digits=1)
end

"""
Calculate the site classification given the fitted model and index age.

# Arguments
- `fitted_model`: The fitted regression model.
- `index_age::Real`: The index age for site classification.

# Returns
- `site_classification::Vector{Real}`: The site classification for each observation.
"""
function site_classification(fitted_model::FittedLinearModel, index_age::Real)
  # Extract the data from the fitted model
  data_age = fitted_model.data
  return site_classification(fitted_model, data_age, index_age)
end

"""
Calculate the dominant height given the site classification, fitted model, and index age.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted regression model.
- `data_age::AbstractDataFrame`: The data frame containing the current age data.
- `index_age::Real`: The index age for the calculation.
- `site::Vector{<:Real}`: The site classification values.

# Returns
- `hdom::Vector{Real}`: The dominant height for each observation.
"""
function hdom_classification(fitted_model::FittedLinearModel, data_age::AbstractDataFrame, index_age::Real, site::Vector{<:Real})
  Δ = _calculate_delta(fitted_model, data_age, index_age)
  (yname, xname, qname...) = propertynames(fitted_model.data)
  site_data = deepcopy(data_age[!, [yname, xname, qname...]])
  site_data[!, yname] .= site
  site_data[!, xname] .= index_age
  # Calculate the site classification
  y = fitted_model.formula.lhs
  site = modelcols(y, site_data)
  hdom = site .- Δ # Invert the delta for hdom classification
  if isa(y, FunctionTerm)
    hdom = _predict(hdom, data_age[!, xname], fitted_model.σ², nameof(y.f))
  end
  return round.(hdom, digits=1)
end

"""
Calculate the site table given a fitted model, index age, and height increment.

# Arguments
- `fitted_model::FittedLinearModel`: The fitted regression model.
- `index_age::Real`: The index age for site classification.
- `hi::Real`: The height increment for site classification.

# Returns
- `SiteAnalysis`: A struct containing site table and site plot.
"""
function site_table(fitted_model::FittedLinearModel, index_age::Real, hi::Real)
  # Extract property names from the fitted model's data
  (hd, age, q...) = propertynames(fitted_model.data[!, :])
  # Calculate site classification
  site = site_classification(fitted_model, index_age)
  # Get unique and sorted site classes
  sites = _class_center.(site, hi) |> unique |> sort
  # Get unique and sorted ages
  ages = fitted_model.data[:, age] |> unique |> sort
  # Repeat ages for each site class
  repeated_ages = repeat(ages, outer=length(sites))
  # Create repeated site classes
  repeated_sites = vcat([fill(s, length(ages)) for s in sites]...)
  # Create initial DataFrame with repeated ages and sites
  site_table = DataFrame(age => repeated_ages, hd => repeated_sites)
  # Predict the dominant heights for the site classification
  hdom_predict = hdom_classification(fitted_model, site_table, index_age, repeated_sites)
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
    group = categorical(repeated_sites, levels = sort(sites, rev=true)),
    xlabel = age,
    ylabel = hd,
    title = "Site Classification Plot",
    legend = :outertopright,
    markerstrokewidth = 0,
    seriesalpha = 0.6,
    framestyle = :box,
    margin = 3mm,
    color_palette = cgrad(:darktest, categorical = true),
    tick_direction = :out,
    grid = :none,
    titlefont = 10,
    fontfamily = "times",
    guidefontsize = 9,
    legendfontsize = 7
  )

  return SiteAnalysis(site_table, site_plot)

end

function site_table(fitted_model::FittedLinearModel, index_age::Real)
  # Calculate site classification
  site = site_classification(fitted_model, index_age)
  # Compute the class breadth
  h = _amplitude(site)
  k = _sturges(length(site))
  hi = _class_breadth(h, k)
  return site_table(fitted_model, index_age, hi)
end