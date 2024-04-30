_sturges(n::Int) = ceil(Integer, log2(n)) + 1

_class_center(x::Real, hi::Real) = round(x / hi) * hi + (hi / 2)

_amplitude(x::Vector) = maximum(x) - minimum(x)

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

_simple_frequency(x::Vector) = map(i -> count(==(i), x), unique(x) |> sort)

function frequency_table(x::Vector{<:Real}, hi::Real) :: DataFrame
  n = length(x)
  cc = _class_center.(x, hi)
  Xi = unique(cc) |> sort
  LI = Xi .- hi / 2
  LS = Xi .+ hi / 2
  fi = _simple_frequency(cc)
  Fi = cumsum(fi)
  fri = (fi ./ n) .* 100
  Fri = cumsum(fri)
  DataFrame(LI = LI, Xi = Xi, LS = LS, fi = fi, Fi = Fi, fri = fri, Fri = Fri)
end

function frequency_table(x::Vector{<:Real})
  h = _amplitude(x)
  k = _sturges(length(x))
  hi = _class_breadth(h, k)
  frequency_table(x, hi)
end

function frequency_table(g::S, x::S, data::AbstractDataFrame) where S <: Symbol
  combine(groupby(data, g)) do df
    frequency_table(df[:, x])
  end
end

function frequency_table(g::S, x::S, hi::Real, data::AbstractDataFrame) where S <: Symbol
  combine(groupby(data, g)) do df
    frequency_table(df[:, x], hi)
  end
end

function diametric_table(x::Vector{<:Real}, hi::Real; plot_area::Real=1.0) 
  ftable = frequency_table(x, hi)
  ftable.g = (π / 40000) .* ftable.Xi.^2
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

function diametric_table(x::Vector{<:Real}; plot_area::Real=1.0)
  h = _amplitude(x)
  k = _sturges(length(x))
  hi = _class_breadth(h, k)
  diametric_table(x, hi, plot_area=plot_area)
end

function diametric_table(g::S, x::S, data::AbstractDataFrame; plot_area::Real=1.0) where S <: Symbol
  combine(groupby(data, g)) do df
    diametric_table(df[:, x], plot_area=plot_area)
  end
end

function diametric_table(g::S, x::S, hi::Real, data::AbstractDataFrame; plot_area::Real=1.0) where S <: Symbol
  combine(groupby(data, g)) do df
    diametric_table(df[:, x], hi, plot_area=plot_area)
  end
end