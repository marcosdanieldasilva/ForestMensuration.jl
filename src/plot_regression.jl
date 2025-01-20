function _evaluate_regression_type(model::LinearModel)
  n = length(model.data)  # Number of columns in data    
  if n == 2
    return :simple
  elseif n >= 3 && model.data[3] isa CategoricalVector
    return :simple
  else
    return :multiple
  end
end

function _rgba_to_string(rgba::ColorTypes.RGBA{Float64})
  r, g, b, a = red(rgba), green(rgba), blue(rgba), 0.6
  return "rgba($r, $g, $b, $a)"
end

function _calculate_qq(x::Vector{<:Real})
  dist = fit_mle(Normal, x)
  qq = qqbuild(dist, x)
  return DataFrame(qx=qq.qx, qy=qq.qy)
end

function _graph_table(model::LinearModel)
  y = model.data[1]
  # Predicted values from the model
  ŷ = predict(model)
  # Residuals: the difference between observed and predicted values
  resid = y - ŷ
  # Fit the residuals to a normal distribution and build the Q-Q plot data
  dist = fit_mle(Normal, resid)
  qq = qqbuild(dist, resid)
  # Create a DataFrame combining the original data with the residuals and predicted values
  data = hcat(DataFrame(model.data), DataFrame(pred=ŷ, resid=resid), makeunique=true)
  # sort the data by rediduals values
  sort!(data, size(data, 2))
  # Insert the Q-Q plot values as new columns, ensuring unique column names
  insertcols!(data, :qx => qq.qx, :qy => qq.qy, makeunique=true)
  sort!(data, 2)
  return data
end

"""
    plot_regression(model::LinearModel)
  
The `plot_regression` function generates four essential diagnostic plots to analyze the performance 
  and validity of a linear regression model in Julia. These plots help in assessing the goodness-of-fit,
   the distribution of residuals, and the assumptions underlying linear regression, such as homoscedasticity
    and normality. The function is particularly useful for ensuring that the regression model meets the 
    assumptions required for reliable inference and prediction.

# Parameters:
- `model::LinearModel`: 
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
function plot_regression(model::LinearModel)
  data = DataFrame(model.data)
  n = size(data, 2)
  rtype = _evaluate_regression_type(model)

  if rtype == :simple && n > 2
    ngroups = groupby(data, 3:n).ngroups
    color_palette = cgrad(:darktest, categorical=true, ngroups)
  elseif rtype == :simple
    color_palette = cgrad(:darktest)[1]
  elseif rtype == :multiple && n > 3
    ngroups = groupby(data, 4:n).ngroups
    color_palette = cgrad(:darktest, categorical=true, ngroups)
  else
    color_palette = cgrad(:darktest)[1]
  end

  data = ForestMensuration._graph_table(model)
  resid = data[:, end-2]
  col_names = propertynames(data)
  x_qqline = [extrema(data[:, end-1])...]
  quantx, quanty = quantile(data[:, end-1], [0.25, 0.75]), quantile(data[:, end], [0.25, 0.75])
  slp = diff(quanty) ./ diff(quantx)
  y_qqline = quanty .+ slp .* (x_qqline .- quantx)
  axis = attr(ticks="outside", showline=true, linewidth=1.5, linecolor="black")

  group_names = rtype == :simple ? col_names[3:n] : col_names[4:n]

  trace3 = PlotlyJS.scatter(
    data,
    x=col_names[end-3],
    y=col_names[end-2],
    group=group_names,
    mode="markers",
    showlegend=false
  )
  trace4 = PlotlyJS.histogram(
    data,
    x=col_names[end-2],
    group=group_names,
    showlegend=false,
    histnorm="probability density",
    nbinsx=length(fit(Histogram, resid).weights)
  )
  trace5 = PlotlyJS.scatter(
    data,
    x=col_names[end-1],
    y=col_names[end],
    group=group_names,
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

  if rtype == :simple

    specs = [Spec(kind="xy") Spec(kind="xy"); Spec(kind="histogram") Spec(kind="xy")]

    fig = make_subplots(
      rows=2, cols=2,
      specs=specs,
      subplot_titles=["Observed vs Fitted Value" "Residuals vs Fitted Value" "Histogram of Residuals" "Normal Q-Q Plot"]
    )

    trace1 = PlotlyJS.scatter(
      data,
      x=col_names[2],
      y=col_names[1],
      group=group_names,
      mode="markers",
      showlegend=n > 2 ? true : false
    )
    trace2 = PlotlyJS.scatter(
      data,
      x=col_names[2],
      y=col_names[end-3],
      group=group_names,
      mode="lines",
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

    map(i -> restyle!(i, marker_color=ForestMensuration._rgba_to_string.(color_palette)), [p1, p2, p3, p4])
    add_trace!(p4, trace6)

    fig = [p1 p2; p3 p4]

    relayout!(
      fig,
      barmode="stack",
      font_size=10,
      font_family="Times New Roman",
      margin=attr(l=40, r=40, t=20, b=40),
      paper_bgcolor="LightSteelBlue"
    )

  else

    specs = [Spec(kind="scene", rowspan=3) Spec(kind="xy"); missing Spec(kind="histogram"); missing Spec(kind="xy")]

    fig = make_subplots(
      rows=3, cols=2,
      column_widths=[0.5, 0.4],
      specs=specs,
      subplot_titles=["Observed vs Fitted Value" "Observed $(col_names[2]) vs Fitted Value" "Observed $(col_names[3]) vs Fitted Value" "Residuals vs Fitted Value" "Histogram of Residuals" "Normal Q-Q Plot"],
      vertical_spacing=0.15,
      horizontal_spacing=0.0
    )

    trace1 = PlotlyJS.scatter3d(
      data,
      x=col_names[3],
      y=col_names[2],
      z=col_names[1],
      group=group_names,
      mode="markers",
      showlegend=n > 3 ? true : false
    )

    trace2 = PlotlyJS.mesh3d(
      data,
      x=col_names[3],
      y=col_names[2],
      z=col_names[end-3],
      group=group_names,
      showlegend=false
    )

    # return fig
    map(t -> add_trace!(fig, t, row=1, col=1), trace1)
    map(t -> add_trace!(fig, t, row=1, col=1), trace2)
    map(t -> add_trace!(fig, t, row=1, col=2), trace3)
    map(t -> add_trace!(fig, t, row=2, col=2), trace4)
    map(t -> add_trace!(fig, t, row=3, col=2), trace5)
    add_trace!(fig, trace6, row=3, col=2)

    restyle!(fig, marker_color=ForestMensuration._rgba_to_string.(color_palette))
    restyle!(fig, color=ForestMensuration._rgba_to_string.(color_palette))

    relayout!(
      fig,
      #p1
      scene=attr(
        zaxis_title=col_names[1],
        yaxis_title=col_names[2],
        xaxis_title=col_names[3],
        axis=axis,
        yaxis=axis,
        xaxis=axis
      ),
      scene_camera=attr(
        eye=attr(x=2.5, y=2.5, z=1.0),
        center=attr(x=0, y=0, z=-0.2)),
      #p2
      xaxis_title="Predicted $(col_names[1])",
      yaxis_title="Residuals",
      xaxis=axis,
      yaxis=axis,
      xaxis2=axis,
      yaxis2=yaxis = attr(ticks="outside", showline=true, linewidth=2, linecolor="black", zeroline=true, zerolinecolor="Grey"),
      xaxis3=axis,
      yaxis3=axis,
      xaxis4=axis,
      yaxis4=axis,
      # p3
      xaxis2_title="Residuals",
      yaxis2_title="Probability Density",
      barmode="stack",
      title_font_size=8,
      font_size=8,
      font_family="Times New Roman",
      margin=attr(l=20, r=10, t=30, b=5),
      paper_bgcolor="LightSteelBlue"
    )
  end

  return fig
end