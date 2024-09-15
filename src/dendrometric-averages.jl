basal_area(d::Real) = (π / 40000) * d^2

function dendrometric_averages(d::Vector; area::Real=1.0) :: DataFrame
  d̅, s = mean_and_std(d)
  d₋, d₊ = d̅ - s, d̅ + s
  g = basal_area.(d)
  dg = √((40000 * mean(g)) / π)
  dw = quantile(d, 0.6)
  dz = √((40000 * median(g)) / π)
  n_tree = round(Int, (100 * area))
  d₁₀₀ = n_tree < length(d) ? mean(partialsort(d, 1:n_tree, rev = true)) : NaN64
  return DataFrame(d₋ = d₋, d̅ = d̅, dg = dg, dw = dw, dz = dz, d₁₀₀ = d₁₀₀, d₊ = d₊)
end

function dendrometric_averages(p::S, d::S, data::AbstractDataFrame; area::Real=1.0) :: DataFrame where S <: Symbol
  return combine(groupby(data, p)) do df
    dendrometric_averages(df[:, d], area = area)
  end
end

function dendrometric_averages(model::TableRegressionModel; area::Real=1.0)
  data = model.mf.data
  var_names = propertynames(data) |> collect
  if length(data) == 2
    result_table = dendrometric_averages(data[2], area = area)
    average_heights = map(i -> prediction(model, DataFrame(var_names[2] => result_table[:, i]))[1], 1:7)
    return hcat(
      result_table,
      DataFrame(
        h₋ = average_heights[1], h̅ = average_heights[2], hg = average_heights[3], 
        hw = average_heights[4], hz = average_heights[5], h₊ = average_heights[6], 
        h₁₀₀ = average_heights[7]
      )
    )
  else
    result_table = combine(groupby(DataFrame(data), var_names[3:end])) do df
      dendrometric_averages(df[:, var_names[2]], area = area)
    end
    s = size(result_table, 2)
    average_heights = hcat(
      map(i -> prediction(model, hcat(result_table[:, 1:s-7], DataFrame(var_names[2] => result_table[:, s-i]))), 6:-1:0)...
    )
    return hcat(result_table, DataFrame(average_heights, [:h₋, :h̅, :hg, :hw, :hz, :h₁₀₀, :h₊]))
  end
end