function _diameter_interpolation(h0::Real, d::Vector{<:Real}, h::Vector{<:Real}) :: Float64
  n = length(h)
  for i in 2:n
      if h0 <= h[i]
          d_prev = d[i - 1]  # Diameter at the lower height
          d_next = d[i]      # Diameter at the higher height
          h_prev = h[i - 1]  # Lower height
          h_next = h[i]      # Higher height
          # Perform linear interpolation
          return d_prev + ((d_next - d_prev) * (h0 - h_prev)) / (h_next - h_prev)
      end
  end
  error("Height h0 is outside the range of heights in the data.")
end

function _cylinder_volume(d::Real, h::Real) :: Float64
  (π / 40000) * d^2 * h
end

function _cone_volume(d::Real, h::Real) :: Float64
  (1 / 3) * (π / 40000) * d^2 * h
end

function _bark_factor(d::Vector{<:Real}, e::Vector{<:Real})
  1 - (sum(e) / sum(d))
end

abstract type CubingMethod end

"""
Smalian Method:
  The Smalian method measures diameters or circumferences at the ends of each section and calculates the total volume by:
  Vt = v0 + Σi=1:n(vi) + vt
  v0 = g0 * l0
  vi = (gi+gi+1)/2 * li
  vt = (1/3) * gn * ln
  Where:
    v0 = volume of the stump;
    vi = volume of intermediate sections;
    vt = volume of the cone;
    g = basal area;
    l = length.
"""
abstract type Smalian <: CubingMethod end

@inline _bole_volume(::Type{<:Smalian}, h::Vector{<:Real}, d::Vector{<:Real}) :: Float64 = map(i -> (π / 40000) * ((d[i]^2 + d[i - 1]^2) / 2) * (h[i] - h[i - 1]), 2:(length(h) - 1)) |> sum

"""
Huber Method:
  The Huber method measures the diameter or circumference at the midpoint of the section, and the volume is determined by:
  v = v0 + Σi=1:n(vi) + vt
  vi = gi * li
  Where:
    v0 = volume of the stump;
    vi = volume of intermediate sections;
    vt = volume of the cone;
    g = basal area;
    l = length.
"""
abstract type Huber <: CubingMethod end

@inline _bole_volume(::Type{<:Huber}, h::Vector{<:Real}, d::Vector{<:Real}) :: Float64 =  _cylinder_volume.(map(i -> d[i], 2:2:(length(d) - 1)), map(i -> h[i] - h[i - 2], 3:2:length(h))) |> sum

"""
Newton Method:
  The Newton method involves measuring at 3 positions along each section (at the ends and in the middle of the logs). Therefore, it is a more laborious method than the others, but the estimated volume will be more accurate.
  v = v0 + Σi=1:n(vi) + vt
  vi = (gi + gm + gi+1)/2 * li
  Where:
    v0 = volume of the stump;
    vi = volume of intermediate sections;
    vt = volume of the cone;
    g = basal area;
    gm = basal area at the midpoint of the section;
    l = length.
"""
abstract type Newton <: CubingMethod end

@inline _bole_volume(::Type{<:Newton}, h::Vector{<:Real}, d::Vector{<:Real}) :: Float64 = map(i -> (π / 240000) * (h[i + 2] - h[i]) * (d[i]^2 + 4 * d[i + 1]^2 + d[i + 2]^2), 1:2:(length(h) - 3)) |> sum

"""
Calculate tree cubage using Smalian, Newton, or Huber methods.

The methods involve dividing the tree trunk into n sections (logs). In each section, diameters and lengths are measured at positions that vary according to the technique employed. Thus, the volume of the sections and the total volume are determined by summing the volume of the sections. Determination can be carried out on felled trees or standing trees using equipment such as the Bitterlich relascope.

- Form Factor (ff):
  The form factor is a reduction factor used to determine the volume of a standing tree. It is the ratio of the rigorous volume of the tree to the volume of a reference cylinder. The volume is calculated by the product of the form factor, basal area, and total height.
  v = g + h + f
  Where:
    g = basal area;
    h = total height;
    f = form factor, either natural or artificial.

- Artificial Form Factor (aff: f₁,₃):
  For the calculation of the artificial form factor, the volume of the reference cylinder will have a diameter equal to the tree's DBH.
  f1,3 = Rigorous Vol / Cylinder Vol 1,3
  Where:
    Rigorous Vol = total volume determined by one of the methods: Smalian, Huber, or Newton;
    Cylinder Vol 1,3 = volume of a cylinder with height and diameter equal to the total height and DBH of the tree.

- Natural Form Factor (nff: f₀,₁ₕ):
  For the calculation of the natural form factor, the volume of the reference cylinder will have a diameter equal to the diameter taken at 1/10 of the total height.
  f0,1h = Rigorous Vol / Cylinder Vol 0,1
  Where:
    Rigorous Vol = total volume determined by one of the methods: Smalian, Huber, or Newton;
    Cylinder Vol 0,1 = volume of a cylinder with height equal to the total height of the tree and diameter taken at 1/10 of the total height.
    Interpolate diameter at a given height using linear interpolation.

- Form Quotient (qf: Q):
    The natural decrease in diameter along the trunk defines the so-called form quotient, which is a ratio between diameters. An example of a form quotient is the Schiffel form quotient, given by:
  
    Q = D(1/2H) / DBH, Where
      Q < 1
      D(1/2H) = diameter measured at half the total height of the tree.
    
    Similar to the form factor, the volume of a tree, with or without bark, can be obtained by multiplying the volume of a cylinder by the average form quotient, suitable for the species and the desired volume to be estimated.
    """
function cubage(method::Type{<:CubingMethod}, h::Vector{<:Real}, d::Vector{<:Real}) :: DataFrame
    # Find position where h = 1.3
    idx = findfirst(isequal(1.3), h)
    
    if idx === nothing
        error("Value of h = 1.3 not found.")
    end
    # diameter at basal height and total height of the tree
    dbh = d[idx]
    ht = h[end]
    # Calculate volumes
    cylinder_volume = _cylinder_volume(d[begin], h[begin])
    bole_volume = _bole_volume(method, h, d)
    cone_volume = _cone_volume(d[end - 1], (ht - h[end - 1]))
    vt = cylinder_volume + bole_volume + cone_volume
    # Calculate the Form Factors
    artificial_form_factor = vt / _cylinder_volume(dbh, ht)
    natural_form_factor = vt / _cylinder_volume(_diameter_interpolation(ht * 0.1, d, h), h[end])
    quotient_form = _diameter_interpolation(ht * 0.5, d, h) / dbh
    # Create DataFrame with results
    DataFrame(
        vt = vt, 
        vc = bole_volume,
        dbh = dbh,
        ht = ht,
        aff = artificial_form_factor,
        nff = natural_form_factor,
        qf = quotient_form
    )
end

function cubage(method::Type{<:CubingMethod}, h::Vector{<:Real}, d::Vector{<:Real}, e::Vector{<:Real})
  cubage_table = cubage(method, h, d)
  # bark factor
  insertcols!(cubage_table, :bf => _bark_factor(d, e))
  # total volume without bark
  insertcols!(cubage_table, :vtwb => cubage_table.vt .* cubage_table.bf.^2)
  # comercial volume without bark
  insertcols!(cubage_table, :vcwb => cubage_table.vc .* cubage_table.bf.^2)
  return cubage_table
end

function cubage(method::Type{<:CubingMethod}, tree::Symbol, h::Symbol, d::Symbol, data::AbstractDataFrame)
    combine(groupby(data, tree)) do df
      cubage(method, df[:, h], df[:, d])
    end
end

function cubage(method::Type{<:CubingMethod}, tree::Symbol, h::Symbol, d::Symbol, e::Symbol, data::AbstractDataFrame)
  combine(groupby(data, tree)) do df
    cubage(method, df[:, h], df[:, d], df[:, e])
  end
end