# unit aliases, kept local for brevity; the underlying types come from ForestCore
const Len = Length   # Length (d, h)
const Vol = Volume   # Volume (v)

# Default units assumed when a measurement is supplied as a plain number.
# Every `Real` method in this package promotes its arguments with these units and
# then forwards to the corresponding `Unitful` method, so unitless and unitful
# workflows always return the same quantities.
const DUNIT = u"cm"    # diameters and bark thicknesses
const HUNIT = u"m"     # heights and log lengths
const VUNIT = u"m^3"   # volumes
const AUNIT = u"ha"    # plot areas

# Promotes a plain number (or vector of numbers) to a quantity in the given unit and
# leaves quantities untouched, so keyword arguments can be normalised uniformly.
_withunit(x::Union{Real,AbstractVector{<:Real}}, u::Units) = x * u
_withunit(x::Union{Quantity,AbstractVector{<:Quantity}}, ::Units) = x
_withunit(::Nothing, ::Units) = nothing
