function _rgba_to_string(rgba::ColorTypes.RGBA{Float64})
  r, g, b, a = red(rgba), green(rgba), blue(rgba), 0.6
  return "rgba($r, $g, $b, $a)"
end

function _calculate_qq(x::Vector{<:Real})
  dist = fit_mle(Normal, x)
  qq = qqbuild(dist, x)
  return DataFrame(qx=qq.qx, qy=qq.qy)
end

function _graph_table(model::TableRegressionModel)
  y = model.mf.data[1]
  # Predicted values from the model
  ŷ = prediction(model)
  # Residuals: the difference between observed and predicted values
  resid = y - ŷ
  # Fit the residuals to a normal distribution and build the Q-Q plot data
  dist = fit_mle(Normal, resid)
  qq = qqbuild(dist, resid)
  # Create a DataFrame combining the original data with the residuals and predicted values
  data = hcat(DataFrame(model.mf.data), DataFrame(pred=ŷ, resid=resid), makeunique=true)
  # sort the data by rediduals values
  sort!(data, size(data, 2))
  # Insert the Q-Q plot values as new columns, ensuring unique column names
  insertcols!(data, :qx => qq.qx, :qy => qq.qy, makeunique=true)
  sort!(data, 2)
  return data
end

"""
  plot_regression(model::RegressionModel)
  
The `plot_regression` function generates four essential diagnostic plots to analyze the performance 
  and validity of a linear regression model in Julia. These plots help in assessing the goodness-of-fit,
   the distribution of residuals, and the assumptions underlying linear regression, such as homoscedasticity
    and normality. The function is particularly useful for ensuring that the regression model meets the 
    assumptions required for reliable inference and prediction.

# Parameters:
- `model::TableRegressionModel`: 
  The fitted linear regression model. This model is analyzed to generate the diagnostic plots.

# Functionality:
- **Observed vs Fitted Values**: 
   This scatter plot compares the observed data points with the values predicted by the model. It is 
   useful for visually assessing the accuracy of the model's predictions.

- **Residuals vs Fitted Values**: 
   This plot shows the residuals (errors) against the fitted values. Ideally, residuals should be randomly
    scattered around zero without any apparent pattern, indicating a good model fit.

- **Histogram of Residuals**: 
   This histogram displays the distribution of residuals. A normal distribution of residuals is one of the
    key assumptions of linear regression, and this plot helps verify it.

- **Normal Q-Q Plot**: 
   This plot compares the quantiles of the residuals to the quantiles of a normal distribution. It is used
    to check the assumption of normality in the residuals.

The function automatically applies color schemes and adjusts plot aesthetics to enhance the clarity of the
   visualizations. The plots are combined into a single layout for easy comparison and interpretation. 
   Additionally, the function adjusts dynamically based on the number of groups in the data, making it 
   versatile for different datasets and models.

# Examples:
```julia
# Generate diagnostic plots for a fitted regression model
plots = plot_regression(model)
```
"""
function plot_regression(model::RegressionModel)
  data = DataFrame(model.mf.data)
  n = size(data, 2)
  if n > 2
    ngroups = groupby(data, 3:n).ngroups
    color_palette = cgrad(:darktest, categorical=true, ngroups)
  else
    color_palette = cgrad(:darktest)[1]
  end
  data = _graph_table(model)
  resid = data[:, end-2]
  col_names = propertynames(data)
  x_qqline = [extrema(data[:, end-1])...]
  quantx, quanty = quantile(data[:, end-1], [0.25, 0.75]), quantile(data[:, end], [0.25, 0.75])
  slp = diff(quanty) ./ diff(quantx)
  y_qqline = quanty .+ slp .* (x_qqline .- quantx)
  axis = attr(ticks="outside", showline=true, linewidth=2, linecolor="black")

  trace1 = PlotlyJS.scatter(
    data,
    x=col_names[2],
    y=col_names[1],
    group=col_names[3:n],
    mode="markers",
    showlegend=n > 2 ? true : false
  )
  trace2 = PlotlyJS.scatter(
    data,
    x=col_names[2],
    y=col_names[end-3],
    group=col_names[3:n],
    mode="lines",
    showlegend=false
  )
  trace3 = PlotlyJS.scatter(
    data,
    x=col_names[end-3],
    y=col_names[end-2],
    group=col_names[3:n],
    mode="markers",
    showlegend=false
  )
  trace4 = PlotlyJS.histogram(
    data,
    x=col_names[end-2],
    group=col_names[3:n],
    showlegend=false,
    histnorm="probability density",
    nbinsx=length(fit(Histogram, resid).weights)
  )
  trace5 = PlotlyJS.scatter(
    data,
    x=col_names[end-1],
    y=col_names[end],
    group=col_names[3:n],
    mode="markers",
    showlegend=false
  )
  trace6 = PlotlyJS.scatter(
    ;
    x=x_qqline,
    y=y_qqline,
    mode="lines",
    marker_color="Grey",
    showlegend=false
  )

  p1 = PlotlyJS.plot(
    trace1,
    Layout(
      xaxis_title=col_names[2],
      yaxis_title=col_names[1],
      title="Observed vs Fitted Value",
      xaxis=axis,
      yaxis=axis
    )
  )

  [add_trace!(p1, t) for t in trace2]

  p2 = PlotlyJS.plot(
    trace3,
    Layout(
      xaxis_title="Predicted $(col_names[1])",
      yaxis_title="Residuals",
      title="Residuals vs Fitted Value",
      size=8,
      xaxis=axis,
      yaxis=attr(ticks="outside", showline=true, linewidth=2, linecolor="black", zeroline=true, zerolinecolor="Grey")
    )
  )

  p3 = PlotlyJS.plot(
    trace4,
    Layout(
      xaxis_title="Residuals",
      yaxis_title="Probability Density",
      title="Histogram of Residuals",
      xaxis=axis,
      yaxis=axis
    )
  )

  p4 = PlotlyJS.plot(
    trace5,
    Layout(
      xaxis_title="Theoretical Quantiles",
      yaxis_title="Empirical Quantiles",
      title="Normal Q-Q Plot",
      xaxis=axis,
      yaxis=axis
    )
  )

  map(i -> restyle!(i, marker_color=_rgba_to_string.(color_palette)), [p1, p2, p3, p4])
  add_trace!(p4, trace6)

  p = [p1 p2; p3 p4]

  relayout!(
    p,
    barmode="stack",
    font_size=10,
    font_family="Times New Roman",
    margin=attr(l=40, r=40, t=20, b=40),
    paper_bgcolor="LightSteelBlue"
  )

  return p
end
