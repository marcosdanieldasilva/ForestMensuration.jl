"""
    smalian(L::Len, dbase::Len, dtop::Len)

Calculates the volume of a single log using Smalian's method.

# Arguments
- `L::Len`: Log length.
- `dbase::Len`: Diameter at the large end of the log (base).
- `dtop::Len`: Diameter at the small end of the log (top).

# Returns
- `Quantity`: The volume of the log in cubic units (e.g., cubic meters `m^3` or cubic feet `ft^3`).

# Mathematical basis
The method assumes the log resembles a frustum of a paraboloid. The volume is calculated 
by averaging the cross-sectional areas of the two ends and multiplying by the length:

```math
v = \\frac{B + b}{2} L

```

Where:

* `B` = cross-sectional area at the large end of the log.
* `b` = cross-sectional area at the small end of the log.
* `L` = log length.

# Examples

```julia-repl
julia> smalian(3.0u"m", 30.0u"cm", 25.0u"cm")
0.17965982987716633 m^3

julia> smalian(10.0u"ft", 12.0u"inch", 10.0u"inch")
6.654067773228381 ft^3

```

"""
function smalian(L::Len, dbase::Len, dtop::Len)
  L <= zero(L) && throw(DomainError(L, "Log length must be greater than zero."))
  B = basalarea(dbase)
  b = basalarea(dtop)
  ((B + b) / 2) * L
end

"""
    smalian(h::AbstractVector{<:Len}, d::AbstractVector{<:Len})

Calculates the total bole volume from a continuous profile of cumulative heights.

# Arguments

* `h::AbstractVector{<:Len}`: Vector of cumulative heights from the base of the tree.
* `d::AbstractVector{<:Len}`: Vector of diameters corresponding to each height.

# Returns

* `Quantity`: The total volume of the bole in cubic units (e.g., cubic meters `m^3` or cubic feet `ft^3`).

# Mathematical basis

The total volume is the sum of the volumes of all individual sections, where each section is
calculated using Smalian's method. The length of each section is the difference between
consecutive cumulative heights:

```math
v = \\sum_{i=2}^{n} \\frac{g_{i-1} + g_i}{2} (h_i - h_{i-1})

```

Where:

* `g_i` = cross-sectional area at cumulative height `h_i`.
* `h_i` = cumulative height at measurement position `i`.
* `n` = total number of measurement points along the bole.

# Examples

```julia-repl
julia> hvec = [0.1, 1.3, 3.3, 5.3] .* u"m";
julia> dvec = [30.0, 25.0, 18.0, 10.0] .* u"cm";
julia> smalian(hvec, dvec)
0.17969909978533616 m^3

```

"""
function smalian(h::AbstractVector{<:Len}, d::AbstractVector{<:Len})
  length(h) == length(d) || throw(DimensionMismatch("Height and diameter vectors must have the same length."))
  sum(i -> smalian(h[i] - h[i-1], d[i-1], d[i]), 2:length(h))
end

"""
    smalian(L::Len, d::AbstractVector{<:Len})

Calculates the total bole volume using Smalian's method with relative (equal-length) sections.

# Arguments
- `L::Len`: Total length of the bole.
- `d::AbstractVector{<:Len}`: Vector of diameters measured at equal intervals along the bole (including base and top).

# Returns
- `Quantity`: The total volume in cubic units (e.g., cubic meters `m^3` or cubic feet `ft^3`).

# Mathematical basis
When a bole is divided into sections of equal length, Smalian's formula can be mathematically factored to improve performance. The constant section length is multiplied by the sum of the intermediate cross-sectional areas and the average of the two end areas:

```math
v = \\frac{L}{n} \\left( \\frac{g_1 + g_{n+1}}{2} + \\sum_{i=2}^{n} g_i \\right)

```

Where:

* `L` = total length of the bole.
* `n` = number of equal-length sections.
* `g_i` = cross-sectional area at measurement position `i`.

# Examples

```julia-repl
julia> L = 25.0u"m"
julia> dvec = [30.0, 28.0, 25.0, 22.0, 19.0, 16.0, 13.0, 10.0, 7.0, 4.0, 0.0] .* u"cm"
julia> smalian(L, dvec)
0.6467753875577987 m^3

```

"""
function smalian(L::Len, d::AbstractVector{<:Len})
  L <= zero(L) && throw(DomainError(L, "Total length must be greater than zero."))
  nsections = length(d) - 1
  nsections < 1 && throw(ArgumentError("Diameter vector must have at least two measurements."))
  lsec = L / nsections
  gends = (basalarea(d[begin]) + basalarea(d[end])) / 2
  gmid = sum(basalarea, @view d[2:(end-1)])
  (gends + gmid) * lsec
end

"""
    hohenadl(L::Len, d::AbstractVector{<:Len})

Calculates the total bole volume using Hohenadl's method based on relative section centers.

# Arguments

* `L::Len`: Total length of the bole.
* `d::AbstractVector{<:Len}`: Vector of diameters measured at the exact center of each equal-length relative section (e.g., 10%, 30%, 50%, 70%, 90% of total length).

# Returns

* `Quantity`: The total volume in cubic units (e.g., cubic meters `m^3` or cubic feet `ft^3`).

# Mathematical basis

Hohenadl's method is a special case of Huber's method applied to sections of equal length. The total volume is the product of the constant section length and the sum of the cross-sectional areas at the midpoint of each section:

```math
v = \\frac{L}{n} \\sum_{i=1}^{n} g_i

```

Where:

* `L` = total length of the bole.
* `n` = number of sections.
* `g_i` = cross-sectional area at the midpoint of section `i`.

# Examples

```julia-repl
julia> L = 25.0u"m"
julia> dvec = [30.0, 28.0, 25.0, 22.0, 19.0, 16.0, 13.0, 10.0, 7.0, 4.0, 0.0] .* u"cm"
julia> hohenadl(L, dvec)
0.6683024372181924 m^3

```

"""
function hohenadl(L::Len, d::AbstractVector{<:Len})
  L <= zero(L) && throw(DomainError(L, "Total length must be greater than zero."))
  nsections = length(d)
  nsections < 1 && throw(ArgumentError("Diameter vector cannot be empty."))
  lsec = L / nsections
  gmid = sum(basalarea, d)
  gmid * lsec
end

"""
    huber(L::Len, dmid::Len)

Calculates the volume of a single log using Huber's method.

# Arguments

* `L::Len`: Log length.
* `dmid::Len`: Diameter at the exact midpoint of the log.

# Returns

* `Quantity`: The volume of the log in cubic units (e.g., cubic meters `m^3` or cubic feet `ft^3`).

# Mathematical basis

The method assumes the log resembles a paraboloid frustum. The volume is calculated
by taking the cross-sectional area at the log's midpoint and multiplying by its length:

```math
v = B_{1/2} L

```

Where:

* `B_{1/2}` = cross-sectional area at the midpoint of the log length.
* `L` = log length.

# Examples

```julia-repl
julia> huber(3.0u"m", 27.5u"cm")
0.1781872083207961 m^3

julia> huber(10.0u"ft", 11.0u"inch")
6.599526234103559 ft^3

```

"""
function huber(L::Len, dmid::Len)
  L <= zero(L) && throw(DomainError(L, "Log length must be greater than zero."))
  bmid = basalarea(dmid)
  bmid * L
end

"""
    newton(L::Len, dbase::Len, dmid::Len, dtop::Len)

Calculates the volume of a single log using Newton's method.

# Arguments

* `L::Len`: Log length.
* `dbase::Len`: Diameter at the large end of the log (base).
* `dmid::Len`: Diameter at the exact midpoint of the log.
* `dtop::Len`: Diameter at the small end of the log (top).

# Returns

* `Quantity`: The volume of the log in cubic units (e.g., cubic meters `m^3` or cubic feet `ft^3`).

# Mathematical basis

Newton's formula (based on Simpson's rule) is the most accurate of the standard
log scaling methods. It calculates volume by weighting the cross-sectional areas:

```math
v = \\frac{B + 4B_{1/2} + b}{6} L

```

Where:

* `B` = cross-sectional area at the large end of the log.
* `B_{1/2}` = cross-sectional area at the midpoint of the log length.
* `b` = cross-sectional area at the small end of the log.
* `L` = log length.

# Examples

```julia-repl
julia> newton(3.0u"m", 30.0u"cm", 27.5u"cm", 25.0u"cm")
0.1786780821729195 m^3

julia> newton(10.0u"ft", 12.0u"inch", 11.0u"inch", 10.0u"inch")
6.617706747145165 ft^3

```

"""
function newton(L::Len, dbase::Len, dmid::Len, dtop::Len)
  L <= zero(L) && throw(DomainError(L, "Log length must be greater than zero."))
  B = basalarea(dbase)
  bmid = basalarea(dmid)
  b = basalarea(dtop)
  ((B + 4 * bmid + b) / 6) * L
end

"""
    cylindervolume(h::Len, d::Len)

Calculates the volume of a cylinder, used to estimate the volume (v0) of the tree stump remaining after clear-cutting.

# Arguments
- `h::Len`: The height of the cylinder.
- `d::Len`: The diameter of the cylinder.

# Returns
- `Quantity`: The volume of the cylinder in cubic units (e.g., cubic meters `m^3` or cubic feet `ft^3`).

# Mathematical basis
The volume of a cylinder is the product of its cross-sectional area and its height:

```math
v = \\frac{\\pi d^2}{4} h

```

# Examples

```julia-repl
julia> cylindervolume(0.15u"m", 30.0u"cm")
0.010602875205865558 m^3

julia> cylindervolume(0.5u"ft", 12.0u"inch")
0.39269908169872414 ft^3

```

"""
function cylindervolume(h::Len, d::Len)
  h <= zero(h) && throw(DomainError(h, "The height must be a positive value."))
  d <= zero(d) && throw(DomainError(d, "The diameter must be a positive value."))
  v = (π / 4) * abs2(d) * h
  unit(d) isa ImperialUnits ? uconvert(u"ft^3", v) : uconvert(u"m^3", v)
end

"""
    conevolume(h::Len, d::Len)

Calculates the volume of a cone, used to estimate the final portion (vn) or tip of the tree.

# Arguments

* `h::Len`: The height of the cone.
* `d::Len`: The diameter at the base of the cone.

# Returns

* `Quantity`: The volume of the cone in cubic units (e.g., cubic meters `m^3` or cubic feet `ft^3`).

# Mathematical basis

The volume of a cone is exactly one-third the volume of a cylinder with the same base area and height:

```math
v = \\frac{\\pi d^2}{12} h

```

# Examples

```julia-repl
julia> conevolume(2.0u"m", 10.0u"cm")
0.005235987755982988 m^3

julia> conevolume(6.0u"ft", 4.0u"inch")
0.17453292519943295 ft^3

```

"""
function conevolume(h::Len, d::Len)
  h <= zero(h) && throw(DomainError(h, "The height must be a positive value."))
  d <= zero(d) && throw(DomainError(d, "The diameter must be a positive value."))
  v = (π / 12) * abs2(d) * h
  unit(d) isa ImperialUnits ? uconvert(u"ft^3", v) : uconvert(u"m^3", v)
end

"""
    diameterinterpolation(h0::Len, h::AbstractVector{<:Len}, d::AbstractVector{<:Len})

Interpolates the diameter at a specific height along the bole using linear interpolation.

# Arguments
- `h0::Len`: The target height where the diameter is to be estimated.
- `h::AbstractVector{<:Len}`: Vector of measured cumulative heights.
- `d::AbstractVector{<:Len}`: Vector of measured diameters corresponding to the heights.

# Returns
- `Quantity`: The interpolated diameter in the same units as `d`.
"""
function diameterinterpolation(h0::Len, h::AbstractVector{<:Len}, d::AbstractVector{<:Len})
  (h0 < h[begin] || h0 > h[end]) && throw(DomainError(h0, "Height is outside the range of measured heights."))
  for i in 2:length(h)
    if h0 <= h[i]
      return d[i-1] + ((d[i] - d[i-1]) * (h0 - h[i-1])) / (h[i] - h[i-1])
    end
  end
end

"""
    heightinterpolation(dlimit::Len, h::AbstractVector{<:Len}, d::AbstractVector{<:Len})

Interpolates the height at a specific commercial diameter limit using linear interpolation.

# Arguments
- `dlimit::Len`: The target commercial diameter limit.
- `h::AbstractVector{<:Len}`: Vector of measured cumulative heights.
- `d::AbstractVector{<:Len}`: Vector of measured diameters.

# Returns
- `Tuple{Quantity, Int}`: A tuple containing the interpolated height and the insertion index.
"""
function heightinterpolation(dlimit::Len, h::AbstractVector{<:Len}, d::AbstractVector{<:Len})
  (dlimit > d[begin] || dlimit < d[end]) && throw(DomainError(dlimit, "Diameter is outside the range of measured diameters."))
  for i in 2:length(d)
    if dlimit <= d[i-1] && dlimit >= d[i]
      hi = h[i-1] + ((h[i] - h[i-1]) * (dlimit - d[i-1])) / (d[i] - d[i-1])
      return (hi, i)
    end
  end
end

"""
    artificialformfactor(vt::Vol, ht::Len, dbh::Len)

Calculates the artificial form factor, representing the ratio of total volume to a reference cylinder based on DBH.

# Arguments
- `vt::Vol`: The total rigorous volume of the tree.
- `ht::Len`: The total height of the tree.
- `dbh::Len`: Diameter at breast height.

# Returns
- `Float64`: The dimensionless artificial form factor.
"""
function artificialformfactor(vt::Vol, ht::Len, dbh::Len)
  vt <= zero(vt) && throw(DomainError(vt, "Volume must be positive."))
  ustrip(uconvert(NoUnits, vt / cylindervolume(ht, dbh)))
end

"""
    naturalformfactor(vt::Vol, ht::Len, h::AbstractVector{<:Len}, d::AbstractVector{<:Len})

Calculates the natural form factor using a reference cylinder based on the diameter at 1/10th of the total height.

# Arguments
- `vt::Vol`: The total rigorous volume of the tree.
- `ht::Len`: The total height of the tree.
- `h::AbstractVector{<:Len}`: Vector of measured heights.
- `d::AbstractVector{<:Len}`: Vector of measured diameters.

# Returns
- `Float64`: The dimensionless natural form factor.
"""
function naturalformfactor(vt::Vol, ht::Len, h::AbstractVector{<:Len}, d::AbstractVector{<:Len})
  vt <= zero(vt) && throw(DomainError(vt, "Volume must be positive."))
  d01 = diameterinterpolation(0.1 * ht, h, d)
  ustrip(uconvert(NoUnits, vt / cylindervolume(ht, d01)))
end

"""
    quotientform(ht::Len, dbh::Len, h::AbstractVector{<:Len}, d::AbstractVector{<:Len})

Calculates the form quotient (e.g., Schiffel form quotient) representing the natural decrease in diameter along the trunk.

# Arguments
- `ht::Len`: The total height of the tree.
- `dbh::Len`: Diameter at breast height.
- `h::AbstractVector{<:Len}`: Vector of measured heights.
- `d::AbstractVector{<:Len}`: Vector of measured diameters.

# Returns
- `Float64`: The dimensionless form quotient.
"""
function quotientform(ht::Len, dbh::Len, h::AbstractVector{<:Len}, d::AbstractVector{<:Len})
  dbh <= zero(dbh) && throw(DomainError(dbh, "Diameter must be positive."))
  d05 = diameterinterpolation(0.5 * ht, h, d)
  ustrip(uconvert(NoUnits, d05 / dbh))
end

"""
    cubage(h::AbstractVector{<:Len}, d::AbstractVector{<:Len}; dlimit::Union{Len,Nothing}=nothing, ht::Union{Len,Nothing}=nothing, dbh::Union{Len,Nothing}=nothing, hdbh::Len=1.3u"m")

Calculates the partitioned volume of a single tree, returning a DataFrame with detailed commercial and residual volumes alongside form factors.

# Arguments
- `h::AbstractVector{<:Len}`: Vector of cumulative heights from the ground.
- `d::AbstractVector{<:Len}`: Vector of diameters corresponding to each height.
- `dlimit::Union{Len,Nothing}`: Commercial diameter limit. If `nothing`, defaults to the penultimate diameter `d[end-1]`.
- `ht::Union{Len,Nothing}`: Total height. If `nothing`, assumes the total height is the last measured position `h[end]`.
- `dbh::Union{Len,Nothing}`: Diameter at Breast Height. If `nothing`, interpolates it automatically at `hdbh`.
- `hdbh::Len`: The standardized height for DBH measurement (defaults to `1.3u"m"`).

# Returns
- `DataFrame`: A single-row DataFrame containing `vt`, `v0`, `vc`, `vr`, `vn`, `d`, `h`, `hc`, `aff`, `nff`, and `qf`.

# Examples
```julia-repl
julia> hvec = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3] .* u"m";
julia> dvec = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9] .* u"cm";
julia> cb = cubage(hvec, dvec)

```

"""
function cubage(h::AbstractVector{<:Len}, d::AbstractVector{<:Len}; dlimit::Union{Len,Nothing}=nothing, ht::Union{Len,Nothing}=nothing, dbh::Union{Len,Nothing}=nothing, hdbh::Len=1.3u"m")
  length(h) == length(d) || throw(DimensionMismatch("Vectors must have the same length."))
  uh = unit(h[begin])
  ud = unit(d[begin])
  hwork = float.(uconvert.(uh, h))
  dwork = float.(uconvert.(ud, d))
  htval = ht === nothing ? hwork[end] : float(uconvert(uh, ht))
  hdbhval = float(uconvert(uh, hdbh))
  if dbh === nothing
    idx = findfirst(x -> x == hdbhval, hwork)
    dbhval = idx === nothing ? diameterinterpolation(hdbhval, hwork, dwork) : dwork[idx]
  else
    dbhval = float(uconvert(ud, dbh))
  end
  dlimval = dlimit === nothing ? dwork[end-1] : float(uconvert(ud, dlimit))
  if ht !== nothing && htval > hwork[end]
    push!(hwork, htval)
    push!(dwork, zero(dwork[begin]))
  end
  length(hwork) >= 3 || throw(ArgumentError("At least 3 measurements are required."))
  v0 = cylindervolume(hwork[begin], dwork[begin])
  if dlimval > maximum(dwork)
    hcidx = length(hwork) - 1
    hc = hwork[hcidx]
    vc = zero(v0)
    vr = smalian(@view(hwork[1:hcidx]), @view(dwork[1:hcidx]))
  else
    if dlimval ∉ dwork && dlimval >= minimum(dwork)
      hi, idx_insert = heightinterpolation(dlimval, hwork, dwork)
      insert!(hwork, idx_insert, hi)
      insert!(dwork, idx_insert, dlimval)
      hcidx = idx_insert
    elseif dlimval ∈ dwork
      hcidx = findlast(x -> x == dlimval, dwork)
    else
      hcidx = length(hwork) - 1
    end
    hc = hwork[hcidx]
    vc = smalian(@view(hwork[1:hcidx]), @view(dwork[1:hcidx]))
    vr = hcidx < length(hwork) - 1 ? smalian(@view(hwork[hcidx:(end-1)]), @view(dwork[hcidx:(end-1)])) : zero(vc)
  end
  vn = conevolume(htval - hwork[end-1], dwork[end-1])
  vt = v0 + vc + vr + vn
  aff = artificialformfactor(vt, htval, dbhval)
  nff = naturalformfactor(vt, htval, hwork, dwork)
  qf = quotientform(htval, dbhval, hwork, dwork)
  DataFrame(vt=vt, v0=v0, vc=vc, vr=vr, vn=vn, d=dbhval, h=htval, hc=hc, aff=aff, nff=nff, qf=qf)
end

"""
    barkfactor(d::AbstractVector{<:Len}, e::AbstractVector{<:Len})

Calculates the bark factor (k), used to estimate the volume without bark.

The bark factor is calculated considering the ratio of double bark thickness
to total diameter. The function automatically handles and converts any combination of length units.

# Arguments

* `d::AbstractVector{<:Len}`: A vector of diameters at breast height. The diameters must be strictly positive values.
* `e::AbstractVector{<:Len}`: A vector of double bark thicknesses. The bark thicknesses must be non-negative values.

# Returns

* `Float64`: The dimensionless bark factor, representing the proportion of the diameter without bark.

# Mathematical basis

The bark factor estimates the ratio of the diameter without bark to the diameter over bark across the sampled stand:

```math
k = 1 - \\frac{\\sum e}{\\sum d}

```

Where:

* `e` = double bark thickness.
* `d` = diameter over bark.

# Examples

```julia-repl
julia> dvec = [30.8, 14.7, 24.8, 20.0, 21.7, 7.9, 15.8, 12.7, 18.6, 17.0] .* u"cm";
julia> evec = [2.2, 1.1, 1.2, 0.8, 1.1, 0.4, 1.1, 1.3, 1.0, 1.0] .* u"cm";
julia> barkfactor(dvec, evec)
0.9391304347826087

```

"""
function barkfactor(d::AbstractVector{<:Len}, e::AbstractVector{<:Len})
  any(x -> x <= zero(x), d) && throw(DomainError("Diameters must be positive."))
  any(x -> x < zero(x), e) && throw(DomainError("Bark thickness must be non-negative."))
  1 - uconvert(NoUnits, sum(e) / sum(d))
end

"""
    barkinterpolation(h0::Len, h::AbstractVector{<:Len}, e::AbstractVector{<:Len})

Interpolates the double bark thickness at a specific height along the bole using linear interpolation.

# Arguments
- `h0::Len`: The target height where the bark thickness is to be estimated.
- `h::AbstractVector{<:Len}`: Vector of measured cumulative heights.
- `e::AbstractVector{<:Len}`: Vector of measured double bark thicknesses.

# Returns
- `Quantity`: The interpolated double bark thickness in the same units as `e`.
"""
function barkinterpolation(h0::Len, h::AbstractVector{<:Len}, e::AbstractVector{<:Len})
  (h0 < h[begin] || h0 > h[end]) && throw(DomainError(h0, "Height is outside the range of measured heights."))
  for i in 2:length(h)
    if h0 <= h[i]
      return e[i-1] + ((e[i] - e[i-1]) * (h0 - h[i-1])) / (h[i] - h[i-1])
    end
  end
end

"""
    cubage(h::AbstractVector{<:Len}, d::AbstractVector{<:Len}, e::AbstractVector{<:Len}; dlimit::Union{Len,Nothing}=nothing, ht::Union{Len,Nothing}=nothing, dbh::Union{Len,Nothing}=nothing, hdbh::Len=1.3u"m")

Calculates the partitioned volume of a single tree Over Bark (OB) and Under Bark (UB), returning a comprehensive DataFrame.

# Arguments
- `h::AbstractVector{<:Len}`: Vector of cumulative heights from the ground.
- `d::AbstractVector{<:Len}`: Vector of diameters over bark.
- `e::AbstractVector{<:Len}`: Vector of double bark thicknesses.
- `dlimit::Union{Len,Nothing}`: Commercial diameter limit (evaluated over bark).
- `ht::Union{Len,Nothing}`: Total height.
- `dbh::Union{Len,Nothing}`: Diameter at Breast Height over bark.
- `hdbh::Len`: Standardized height for DBH (defaults to `1.3u"m"`).

# Returns
- `DataFrame`: A single-row DataFrame containing all OB (Over Bark) and UB (Under Bark) volumes, total bark volume (`vbark`), bark factor (`k`), and form factors.

# Examples
```julia-repl
julia> hvec = [0.3, 1.3, 3.3, 5.3] .* u"m";
julia> dvec = [9.0, 7.0, 5.8, 5.1] .* u"cm";
julia> evec = [1.2, 0.8, 0.6, 0.4] .* u"cm";
julia> cb = cubage(hvec, dvec, evec, ht=7.0u"m")

```

"""
function cubage(h::AbstractVector{<:Len}, d::AbstractVector{<:Len}, e::AbstractVector{<:Len}; dlimit::Union{Len,Nothing}=nothing, ht::Union{Len,Nothing}=nothing, dbh::Union{Len,Nothing}=nothing, hdbh::Len=1.3u"m")
  length(h) == length(d) == length(e) || throw(DimensionMismatch("Vectors h, d, and e must have the same length."))
  uh = unit(h[begin])
  ud = unit(d[begin])
  hwork = float.(uconvert.(uh, h))
  dwork = float.(uconvert.(ud, d))
  ework = float.(uconvert.(ud, e))
  htval = ht === nothing ? hwork[end] : float(uconvert(uh, ht))
  hdbhval = float(uconvert(uh, hdbh))
  if dbh === nothing
    idx = findfirst(x -> x == hdbhval, hwork)
    dbhval = idx === nothing ? diameterinterpolation(hdbhval, hwork, dwork) : dwork[idx]
  else
    dbhval = float(uconvert(ud, dbh))
  end
  dlimval = dlimit === nothing ? dwork[end-1] : float(uconvert(ud, dlimit))
  if ht !== nothing && htval > hwork[end]
    push!(hwork, htval)
    push!(dwork, zero(dwork[begin]))
    push!(ework, zero(ework[begin]))
  end
  length(hwork) >= 3 || throw(ArgumentError("At least 3 measurements are required."))
  horig = copy(hwork)
  dorig = copy(dwork)
  eorig = copy(ework)
  if dlimval > maximum(dorig)
    hcidx = length(hwork) - 1
    isallresidual = true
  else
    isallresidual = false
    if dlimval ∉ dwork && dlimval >= minimum(dwork)
      hi, idxinsert = heightinterpolation(dlimval, horig, dorig)
      ei = barkinterpolation(hi, horig, eorig)
      insert!(hwork, idxinsert, hi)
      insert!(dwork, idxinsert, dlimval)
      insert!(ework, idxinsert, ei)
      hcidx = idxinsert
    elseif dlimval ∈ dwork
      hcidx = findlast(x -> x == dlimval, dwork)
    else
      hcidx = length(hwork) - 1
    end
  end
  hc = hwork[hcidx]
  dubwork = dwork .- ework
  any(x -> x < zero(x), dubwork) && throw(DomainError("Double bark thickness exceeds diameter at one or more points."))
  v0ob = cylindervolume(hwork[begin], dwork[begin])
  v0ub = cylindervolume(hwork[begin], dubwork[begin])
  vnob = conevolume(htval - hwork[end-1], dwork[end-1])
  vnub = conevolume(htval - hwork[end-1], dubwork[end-1])
  if isallresidual
    vcob = zero(v0ob)
    vrob = smalian(@view(hwork[1:hcidx]), @view(dwork[1:hcidx]))
    vcub = zero(v0ub)
    vrub = smalian(@view(hwork[1:hcidx]), @view(dubwork[1:hcidx]))
  else
    vcob = smalian(@view(hwork[1:hcidx]), @view(dwork[1:hcidx]))
    vrob = hcidx < length(hwork) - 1 ? smalian(@view(hwork[hcidx:(end-1)]), @view(dwork[hcidx:(end-1)])) : zero(vcob)
    vcub = smalian(@view(hwork[1:hcidx]), @view(dubwork[1:hcidx]))
    vrub = hcidx < length(hwork) - 1 ? smalian(@view(hwork[hcidx:(end-1)]), @view(dubwork[hcidx:(end-1)])) : zero(vcub)
  end
  vtob = v0ob + vcob + vrob + vnob
  vtub = v0ub + vcub + vrub + vnub
  vbark = vtob - vtub
  mask = dorig .> zero(dorig[begin])
  k = barkfactor(@view(dorig[mask]), @view(eorig[mask]))
  affob = artificialformfactor(vtob, htval, dbhval)
  DataFrame(vtob=vtob, vtub=vtub, vbark=vbark, k=k, v0ob=v0ob, v0ub=v0ub, vcob=vcob, vcub=vcub, vrob=vrob, vrub=vrub, vnob=vnob, vnub=vnub, d=dbhval, h=htval, hc=hc, aff=affob)
end

"""
    cubage(id::AbstractVector, h::AbstractVector{<:Len}, d::AbstractVector{<:Len}; kwargs...)

Aggregates and calculates partitioned volumes for multiple trees based on a tree identifier.

# Arguments
- `id::AbstractVector`: Vector of tree identifiers (e.g., integers or strings).
- `h::AbstractVector{<:Len}`: Vector of cumulative heights.
- `d::AbstractVector{<:Len}`: Vector of diameters.
- `kwargs...`: Additional keyword arguments (e.g., `dlimit`, `ht`, `dbh`). These can be provided as single scalar values applied to all trees, or as vectors of the same length as `id` containing repeated tree-level attributes.

# Returns
- `DataFrame`: A consolidated DataFrame with the tree `id` as the first column.
"""
function cubage(id::AbstractVector, h::AbstractVector{<:Len}, d::AbstractVector{<:Len}; kwargs...)
  length(id) == length(h) == length(d) || throw(DimensionMismatch("Vectors id, h, and d must have the same length."))
  ntotal = length(id)
  unique_ids = unique(id)
  results = DataFrame()

  for tree_id in unique_ids
    idx = findall(x -> x == tree_id, id)

    # Dynamically extract tree-specific arguments
    treeargs = Dict{Symbol,Any}()
    for (k, v) in pairs(kwargs)
      if v isa AbstractVector
        length(v) == ntotal || throw(DimensionMismatch("Vector keyword argument \$k must have the same length as the main data vectors."))
        treeargs[k] = v[idx[begin]] # Extracts the repeated scalar value for this specific tree
      else
        treeargs[k] = v # Applies the scalar directly
      end
    end

    tree_df = cubage(@view(h[idx]), @view(d[idx]); treeargs...)
    insertcols!(tree_df, 1, :id => tree_id)
    append!(results, tree_df)
  end

  return results
end

"""
    cubage(id::AbstractVector, h::AbstractVector{<:Len}, d::AbstractVector{<:Len}, e::AbstractVector{<:Len}; kwargs...)

Aggregates and calculates partitioned volumes (Over Bark and Under Bark) for multiple trees based on a tree identifier.

# Arguments
- `id::AbstractVector`: Vector of tree identifiers.
- `h::AbstractVector{<:Len}`: Vector of cumulative heights.
- `d::AbstractVector{<:Len}`: Vector of diameters over bark.
- `e::AbstractVector{<:Len}`: Vector of double bark thicknesses.
- `kwargs...`: Additional keyword arguments, accepted as scalars or full-length vectors.

# Returns
- `DataFrame`: A consolidated DataFrame with the tree `id` and all calculated OB/UB variables.
"""
function cubage(id::AbstractVector, h::AbstractVector{<:Len}, d::AbstractVector{<:Len}, e::AbstractVector{<:Len}; kwargs...)
  length(id) == length(h) == length(d) == length(e) || throw(DimensionMismatch("Vectors id, h, d, and e must have the same length."))
  ntotal = length(id)
  unique_ids = unique(id)
  results = DataFrame()

  for tree_id in unique_ids
    idx = findall(x -> x == tree_id, id)

    treeargs = Dict{Symbol,Any}()
    for (k, v) in pairs(kwargs)
      if v isa AbstractVector
        length(v) == ntotal || throw(DimensionMismatch("Vector keyword argument \$k must have the same length as the main data vectors."))
        treeargs[k] = v[idx[begin]]
      else
        treeargs[k] = v
      end
    end

    tree_df = cubage(@view(h[idx]), @view(d[idx]), @view(e[idx]); treeargs...)
    insertcols!(tree_df, 1, :id => tree_id)
    append!(results, tree_df)
  end

  return results
end

"""
    removeunits(df::AbstractDataFrame; renamecols::Bool=true)

Converts a DataFrame containing Unitful quantities into a plain numeric DataFrame.
Safely ignores dimensionless data (NoUnits) preventing empty parentheses.

# Arguments
- `df::AbstractDataFrame`: The input DataFrame containing columns with physical units.
- `renamecols::Bool`: If `true` (default), appends the unit symbol to the column header.

# Returns
- `DataFrame`: A new DataFrame with plain numeric types, ready for export.
"""
function removeunits(df::AbstractDataFrame; renamecols::Bool=true)
  dfclean = copy(df)

  for colname in names(dfclean)
    col = dfclean[!, colname]
    T = nonmissingtype(eltype(col))

    if T <: Quantity
      u = unit(T)
      # Skip renaming if the column is dimensionless (NoUnits)
      if renamecols && u != NoUnits
        newname = "$(colname) ($(u))"
        rename!(dfclean, colname => newname)
        dfclean[!, newname] = ustrip.(col)
      else
        dfclean[!, colname] = ustrip.(col)
      end
    end
  end

  return dfclean
end

"""
    restoreunits(df::AbstractDataFrame)

Restores Unitful quantities to a plain DataFrame by parsing unit strings located in the column headers.
Safely evaluates mathematical expressions (e.g., m^3) inside the Unitful module scope.

# Arguments
- `df::AbstractDataFrame`: The input DataFrame, typically loaded from a CSV file.

# Returns
- `DataFrame`: A DataFrame with restored Unitful quantities.
"""
function restoreunits(df::AbstractDataFrame)
  dfrestored = copy(df)

  # Regular expression to capture the base name and the unit inside parentheses
  regex = r"^(.*?)\s*\((.*?)\)$"

  for colname in names(dfrestored)
    m = match(regex, colname)

    if m !== nothing
      basename = strip(m.captures[1])
      unitstr = strip(m.captures[2])

      # Handle cases where the parentheses were inadvertently empty
      if unitstr == ""
        rename!(dfrestored, colname => basename)
        continue
      end

      try
        # uparse can return an expression (e.g., for "m^3"). 
        # Core.eval evaluates this expression safely inside the Unitful context.
        u = Core.eval(Unitful, Unitful.uparse(unitstr))

        rename!(dfrestored, colname => basename)
        dfrestored[!, basename] = dfrestored[!, basename] .* u
      catch e
        # Formats the warning outside the macro to prevent scope interpolation bugs
        msg = "Could not parse unit '$(unitstr)' in column '$(colname)'. Error: $e"
        @warn msg
      end
    end
  end

  return dfrestored
end
