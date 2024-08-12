function _diameter_interpolation(h0::Real, h::Vector{<:Real}, d::Vector{<:Real}) :: Float64
  n = length(h)
  for i in 2:n
    if h0 <= h[i]
      d_prev = d[i - 1]  # Diameter at the lower height
      d_next = d[i]    # Diameter at the higher height
      h_prev = h[i - 1]  # Lower height
      h_next = h[i]    # Higher height
      # Perform linear interpolation
      return d_prev + ((d_next - d_prev) * (h0 - h_prev)) / (h_next - h_prev)
    end
  end
  error("Height h0 is outside the range of heights in the data.")
end

# Function to interpolate the height at a given diameter
function _height_interpolation(d_limit::Real, h::Vector{<:Real}, d::Vector{<:Real}) :: Tuple{Float64, Int}
  n = length(d)
  # Find the position where d_limit should be interpolated
  for i in 2:n
    if d_limit <= d[i - 1] && d_limit >= d[i]
      h_prev = h[i - 1]  # Height at the lower diameter
      h_next = h[i]  # Height at the higher diameter
      d_prev = d[i - 1]  # Lower diameter
      d_next = d[i]  # Higher diameter
      # Perform linear interpolation
      interpolated_height = h_prev + ((h_next - h_prev) * (d_limit - d_prev)) / (d_next - d_prev)
      return (interpolated_height, i)
    end
  end

  error("Diameter d_limit is outside the range of diameters in the data.")
end

"""
Calculates the volume of a cylinder, used to estimate the volume (v0) of the tree stump remaining after clear-cutting.

# Arguments
- `h::Real`: The height of the cylinder in meters.
- `d::Real`: The diameter of the cylinder in centimeters.

# Returns
- `Float64`: The volume of the cylinder in cubic meters.
"""
@inline cylinder_volume(h::Real, d::Real) :: Float64 = (π / 40000) * d^2 * h

"""
Calculates the volume of a cone, used to estimate the final portion (vn) of the tree, typically considered to have a conical shape.

# Arguments
- `h::Real`: The height of the cone in meters.
- `d::Real`: The diameter of the cone in centimeters.

# Returns
- `Float64`: The volume of the cone in cubic meters.
"""
@inline cone_volume(h::Real, d::Real) :: Float64 = (1 / 3) * (π / 40000) * d^2 * h

"""
Calculates the bark factor, used to estimate the volume without bark.

The bark factor is used to estimate the volume without bark by considering the ratio of bark thickness to total diameter.

# Arguments
- `d::Vector{<:Real}`: Vector of diameters.
- `e::Vector{<:Real}`: Vector of bark thicknesses in millimeters.

# Returns
- `Float64`: The bark factor.

# Example
# Calculate bark factor
k = bark_factor(d, e)

# Calculate total volume without bark
total_volume_without_bark = total_volume * k^2
"""
@inline bark_factor(d::Vector{<:Real}, e::Vector{<:Real}) :: Float64 = 1 - (sum(e) / sum(d))

"""
Calculate tree bole volume cubic meters using Smalian, Newton, or Huber methods.
The methods involve dividing the tree trunk into n sections (logs). 
In each section, diameters and lengths are measured at positions that vary according to the technique employed.

# Arguments
- `method::Type{<:CubingMethod}`: The method used for cubing (Smalian, Huber, or Newton).
- `h::Vector{<:Real}`: Vector of heights.
- `d::Vector{<:Real}`: Vector of diameters.

# Returns
- `Float64`: The volume of the bole in cubic meters.
"""
@inline bole_volume(method::Type{<:CubingMethod}, h::Vector{<:Real}, d::Vector{<:Real}) :: Float64 = bole_volume(method, h, d)

# Smalian Method for Bole Volume Calculation
@inline function bole_volume(::Type{<:Smalian}, h::Vector{<:Real}, d::Vector{<:Real}, start_idx::Int, end_idx::Int) :: Float64
  return map(i -> (π / 40000) * ((d[i]^2 + d[i - 1]^2) / 2) * (h[i] - h[i - 1]), start_idx:end_idx) |> sum
end
# Huber Method for Bole Volume Calculation
@inline function bole_volume(::Type{<:Huber}, h::Vector{<:Real}, d::Vector{<:Real}, start_idx::Int, end_idx::Int) :: Float64
  cylinder_volume.(map(i -> h[i] - h[i - 2], start_idx + 1:2:end_idx+1), map(i -> d[i], start_idx:2:end_idx)) |> sum
end
# Newton Method for Bole Volume Calculation
@inline function bole_volume(::Type{<:Newton}, h::Vector{<:Real}, d::Vector{<:Real}, start_idx::Int, end_idx::Int) :: Float64
  map(i -> (π / 240000) * (h[i + 2] - h[i]) * (d[i]^2 + 4 * d[i + 1]^2 + d[i + 2]^2), start_idx:2:(end_idx - 2)) |> sum
end

"""
Artificial Form Factor (aff):
For the calculation of the artificial form factor, the volume of the reference cylinder will have a diameter equal to the tree's DBH.

- aff = Rigorous Vol / Cylinder Vol 1.3
Where:
- Rigorous Vol = total volume determined by one of the methods: Smalian, Huber, or Newton;
- Cylinder Vol 1.3 = volume of a cylinder with height and diameter equal to the total height and DBH of the tree.

# Arguments
- `vt::Real`: The total volume of the tree.
- `ht::Real`: The total height of the tree.
- `dbh::Real`: The diameter at breast height of the tree.

# Returns
- `Float64`: The artificial form factor.
"""
@inline artificial_form_factor(vt::Real, ht::Real, dbh::Real) :: Float64 = vt / cylinder_volume(ht, dbh)

"""
Natural Form Factor (nff):
For the calculation of the natural form factor, the volume of the reference cylinder will have a diameter equal to the diameter taken at 1/10 of the total height.

- f0.1h = Rigorous Vol / Cylinder Vol 0.1
Where:
- Rigorous Vol = total volume determined by one of the methods: Smalian, Huber, or Newton;
- Cylinder Vol 0.1 = volume of a cylinder with height equal to the total height of the tree and diameter taken at 1/10 of the total height.
Interpolate diameter at a given height using linear interpolation.

# Arguments
- `vt::Real`: The total volume of the tree.
- `ht::Real`: The total height of the tree.
- `h::Vector{<:Real}`: Vector of heights.
- `d::Vector{<:Real}`: Vector of diameters.

# Returns
- `Float64`: The natural form factor.
"""
@inline natural_form_factor(vt::Real, ht::Real, h::Vector{<:Real}, d::Vector{<:Real}) :: Float64 = vt / cylinder_volume(h[end], _diameter_interpolation(ht * 0.1, h, d))

"""
Form Quotient (qf):
The natural decrease in diameter along the trunk defines the so-called form quotient, which is a ratio between diameters. An example of a form quotient is the Schiffel form quotient, given by:

- Q = D(1/2H) / DBH
Where:
- Q < 1
- D(1/2H) = diameter measured at half the total height of the tree.
  
Similar to the form factor, the volume of a tree, with or without bark, can be obtained by multiplying the volume of a cylinder by the average form quotient, suitable for the species and the desired volume to be estimated.

# Arguments
- `ht::Real`: The total height of the tree.
- `dbh::Real`: The diameter at breast height of the tree.
- `h::Vector{<:Real}`: Vector of heights.
- `d::Vector{<:Real}`: Vector of diameters.

# Returns
- `Float64`: The form quotient.
"""
@inline quotient_form(ht::Real, dbh::Real, h::Vector{<:Real}, d::Vector{<:Real}) :: Float64 = _diameter_interpolation(ht * 0.5, h, d) / dbh

"""
Calculate tree cubage using Smalian, Newton, or Huber methods.
The methods involve dividing the tree trunk into n sections (logs). 
In each section, diameters and lengths are measured at positions that vary according to the technique employed. 
Thus, the volume of the sections and the total volume are determined by summing the volume of the sections. 
Determination can be carried out on felled trees or standing trees using equipment such as the Bitterlich relascope.

# Arguments
- `method::Type{<:CubingMethod}`: The method used for cubing (Smalian, Huber, or Newton).
- `h::Vector{<:Real}`: Vector of heights.
- `d::Vector{<:Real}`: Vector of diameters.
- `d_limit::Union{Float64, Nothing}`: Comercial diameter limit to be used in calculations (optional).
- `dbh::Float64`: Diameter at breast height (default is 1.3 meters).

# Returns
- `DataFrame`: A DataFrame with the calculated volumes and form factors.
"""
function cubage(method::Type{<:CubingMethod}, h::Vector{<:Real}, d::Vector{<:Real}, d_limit::Union{Real, Nothing}=nothing; dbh::Real=1.3) :: DataFrame
  # Find position where h = 1.3
  idx = findfirst(isequal(dbh), h)
  
  if idx === nothing
    error("Value of dbh when h = $dbh not found.")
  end
  
  # total height of the tree
  ht = h[end]

  # If commercial diameter is not provided, set it to the last value in the diameter vector
  if d_limit === nothing
    d_limit = d[end - 1]
  end
  
  if d_limit > maximum(d)
    # If no commercial limit or the limit is greater than any diameter, consider all as residual
    hc = h[end - 1]
    hc_idx = length(h) - 1
    vc = 0.0  # No commercial bole volume
    vr = bole_volume(method, h, d, 2, hc_idx)  # All bole volume is residual
  else
    if d_limit ∉ d && d_limit >= minimum(d)
      # Interpolate to find the corresponding height for the commercial diameter limit
      hi_at_d_limit, hc_idx = _height_interpolation(d_limit, h, d)
      h = deepcopy(h)
      d = deepcopy(d)
      insert!(h, hc_idx, hi_at_d_limit)
      insert!(d, hc_idx, d_limit)
      hc = h[hc_idx]
    elseif  d_limit ∈ d
      hc_idx = findlast(isequal(d_limit), d)
      hc = h[hc_idx]
    else
      hc_idx = length(h) - 1
      hc = h[hc_idx]
    end
    vc = bole_volume(method, h, d, 2, hc_idx)  # Commercial bole volume up to d_limit
    vr = bole_volume(method, h, d, hc_idx + 1, length(h) - 1)  # Residual bole volume above d_limit
  end
  # Calculate other volumes
  v0 = cylinder_volume(h[begin], d[begin])  # Volume of the cylinder at the base
  vn = cone_volume(ht - h[end - 1], d[end - 1])  # Volume of the cone above h[end - 1]
  vt = v0 + vc + vr + vn  # Total volume
  # Calculate the Form Factors
  aff = artificial_form_factor(vt, ht, dbh)
  nff = natural_form_factor(vt, ht, h, d)
  qf = quotient_form(ht, dbh, h, d)
  # Create DataFrame with results
  DataFrame(
    vt = vt,
    v0 = v0,
    vc = vc,
    vr = vr,
    vn = vn,
    dbh = d[idx],
    ht = ht,
    hc = hc,
    aff = aff,
    nff = nff,
    qf = qf
  )
end

# function cubage(method::Type{<:CubingMethod}, h::Vector{<:Real}, d::Vector{<:Real}; dbh::Float64=1.3) :: DataFrame
#   # Find position where h = 1.3
#   idx = findfirst(isequal(dbh), h)
  
#   if idx === nothing
#     error("Value of dbh when h = $dbh not found.")
#   end
#   # diameter at basal height and total height of the tree
#   dbh = d[idx]
#   ht = h[end]
#   hc = h[end - 1]
#   # Calculate volumes
#   v0 = cylinder_volume(h[begin], d[begin])
#   vc = bole_volume(method, h, d)
#   vn = cone_volume(ht - hc, d[end - 1])
#   vt = v0 + vc + vn
#   # Calculate the Form Factors
#   aff = artificial_form_factor(vt, ht, dbh)
#   nff = natural_form_factor(vt, ht, h, d)
#   qf = quotient_form(ht, dbh, h, d)
#   # Create DataFrame with results
#   DataFrame(
#     vt = vt,
#     v0 = v0,
#     vc = vc,
#     vn = vn,
#     dbh = dbh,
#     ht = ht,
#     hc = hc,
#     aff = aff,
#     nff = nff,
#     qf = qf
#   )
# end

"""
Calculate tree cubage including bark factor.
The methods involve dividing the tree trunk into n sections (logs). 
In each section, diameters and lengths are measured at positions that vary according to the technique employed.

# Arguments
- `method::Type{<:CubingMethod}`: The method used for cubing (Smalian, Huber, or Newton).
- `h::Vector{<:Real}`: Vector of heights.
- `d::Vector{<:Real}`: Vector of diameters.
- `e::Vector{<:Real}`: Vector of bark thicknesses.
- `d_limit::Union{Float64, Nothing}`: Comercial diameter limit to be used in calculations (optional).
- `dbh::Float64`: Diameter at breast height (default is 1.3 meters).

# Returns
- `DataFrame`: A DataFrame with the calculated volumes, form factors, and bark-adjusted volumes.
"""
function cubage(method::Type{<:CubingMethod}, h::Vector{<:Real}, d::Vector{<:Real}, e::Vector{<:Real}, d_limit::Union{Real, Nothing}=nothing; dbh::Float64=1.3)
  cubage_table = cubage(method, h, d, d_limit, dbh = dbh)
  # bark factor
  insertcols!(cubage_table, :k => bark_factor(d, e))
  # total volume without bark
  insertcols!(cubage_table, :vtwb => cubage_table.vt .* cubage_table.k.^2)
  # cylinder volume without bark
  insertcols!(cubage_table, :v0wb => cubage_table.v0 .* cubage_table.k.^2)
  # bole volume without bark
  insertcols!(cubage_table, :vcwb => cubage_table.vc .* cubage_table.k.^2)
  # residual volume without bark
  insertcols!(cubage_table, :vrwb => cubage_table.vr .* cubage_table.k.^2)
  # cone volume without bark
  insertcols!(cubage_table, :vnwb => cubage_table.vn .* cubage_table.k.^2)
  return cubage_table
end

"""
Calculate tree cubage using grouped data from a DataFrame.
The methods involve dividing the tree trunk into n sections (logs). 
In each section, diameters and lengths are measured at positions that vary according to the technique employed.

# Arguments
- `method::Type{<:CubingMethod}`: The method used for cubing (Smalian, Huber, or Newton).
- `tree::Symbol`: The symbol representing the tree identifier.
- `h::Symbol`: The symbol representing the heights in the DataFrame.
- `d::Symbol`: The symbol representing the diameters in the DataFrame.
- `data::AbstractDataFrame`: The DataFrame containing the tree data.
- `d_limit::Union{Float64, Nothing}`: Comercial diameter limit to be used in calculations (optional).
- `dbh::Float64`: Diameter at breast height (default is 1.3 meters).

# Returns
- `DataFrame`: A DataFrame with the calculated volumes and form factors for each tree.
"""
function cubage(method::Type{<:CubingMethod}, tree::Symbol, h::Symbol, d::Symbol, data::AbstractDataFrame, d_limit::Union{Real, Nothing}=nothing; dbh::Float64=1.3)
  combine(groupby(data, tree)) do df
    cubage(method, df[:, h], df[:, d], d_limit, dbh = dbh)
  end
end

"""
Calculate tree cubage including bark factor using grouped data from a DataFrame.
The methods involve dividing the tree trunk into n sections (logs). 
In each section, diameters and lengths are measured at positions that vary according to the technique employed.

# Arguments
- `method::Type{<:CubingMethod}`: The method used for cubing (Smalian, Huber, or Newton).
- `tree::Symbol`: The symbol representing the tree identifier.
- `h::Symbol`: The symbol representing the heights in the DataFrame.
- `d::Symbol`: The symbol representing the diameters in the DataFrame.
- `e::Symbol`: The symbol representing the bark thicknesses in the DataFrame.
- `data::AbstractDataFrame`: The DataFrame containing the tree data.
- `d_limit::Union{Float64, Nothing}`: Comercial diameter limit to be used in calculations (optional).
- `dbh::Float64`: Diameter at breast height (default is 1.3 meters).

# Returns
- `DataFrame`: A DataFrame with the calculated volumes, form factors, and bark-adjusted volumes for each tree.
"""
function cubage(method::Type{<:CubingMethod}, tree::Symbol, h::Symbol, d::Symbol, e::Symbol, data::AbstractDataFrame, d_limit::Union{Real, Nothing}=nothing; dbh::Float64=1.3)
  combine(groupby(data, tree)) do df
  cubage(method, df[:, h], df[:, d], df[:, e], d_limit, dbh = dbh)
  end
end