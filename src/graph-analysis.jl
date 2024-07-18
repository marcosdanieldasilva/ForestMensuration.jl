# function _confidence_interval(fitted_model::FittedLinearModel)
#   model_matrix = modelmatrix(fitted_model)
#   R = fitted_model.chol.U
#   variance = deviance(fitted_model) / dof_residual(fitted_model)
#   residvar = ones(size(model_matrix, 2)) * variance
#   tvalue = -quantile(TDist(dof_residual(fitted_model)), 0.025)
#   confidence = tvalue * .√((model_matrix / R).^2 * residvar)
#   return  confidence
# end

using StatsPlots
using RecipesBase
using Distributions
using RecipesPipeline
using Plots.PlotMeasures

@recipe function f(model::FittedLinearModel)
  pred = predict(model)
  resid = residuals(model)
  confidence = confidence_interval(model)
  y_name, x_name = names(model.data)[1:2]
  n = size(model.data, 1)
  grid = [1 / (n - 1):1 / (n - 1):1 - (1 / (n - 1));]
  theoritical_quantiles = quantile.(fit_mle(Normal, resid), grid)
  empirical_quantiles = quantile(resid, grid)
  qqline = [first(theoritical_quantiles), last(theoritical_quantiles)]
  grid --> :none
  tick_direction --> :out
  framestyle --> :box
  markerstrokewidth --> 0
  seriesalpha --> 0.6
  legend --> :outertopright
  (y, x, q...) = eachcol(model.data)
  color_palette --> cgrad(:darktest, categorical = true)
  titlefontize --> 10
  guidefontsize --> 9
  legendfontsize --> 7
  fontfamily --> :times
  size --> (1000, 650)
  margin --> 3mm
  layout := (2, 2)
  if size(model.data, 2) == 2
    y, x = eachcol(model.data)
    @series begin
      seriestype := :scatter
      subplot := 1
      label --> "Observed Values"
      x, y
    end
    @series begin
      ribbon --> confidence
      xaxis --> x_name
      yaxis --> y_name
      title --> "Linear Regression"
      label --> "Fitted Values"
      seriestype := :line
      seriesalpha := 0.7
      line --> 2
      subplot := 1
      x, pred
    end
    @series begin
      seriestype := :scatter
      subplot := 2
      label --> "Observed Values"
      pred, resid
    end
    @series begin
      seriestype := :histogram
      xaxis --> "Residuals"
      yaxis --> "Frequency"
      title --> "Histogram of Residuals"
      label --> "Observed Values"
      linewidth := 0.7
      normalize := :pdf
      ylims := (0, Inf)
      subplot := 3
      resid
    end
    @series begin
      seriestype := :scatter
      label --> "Observed Values"
      subplot := 4
      theoritical_quantiles, empirical_quantiles
    end
  else
    gb = RecipesPipeline._extract_group_attributes(tuple(q...))
    labels, idxs = getfield(gb, 1), getfield(gb, 2)
    ngroups = length(labels)
    for i = 1:ngroups
      color = cgrad(:darktest, categorical = true)[i]
      @series begin
        seriestype := :scatter
        subplot := 1
        label := labels[i]
        color := color
        xvals = x[idxs[i]]
        yvals = y[idxs[i]]
        xvals, yvals
      end
      @series begin
        ribbon --> confidence[idxs[i]]
        seriestype := :line
        seriesalpha := 0.7
        line --> 2
        title --> "Linear Regression"
        subplot := 1
        label := :none
        color := color
        xvals = x[idxs[i]]
        predvals = pred[idxs[i]]
        xvals, predvals
      end
      @series begin
        seriestype := :scatter
        subplot := 2
        label := labels[i]
        color := color
        pred[idxs[i]], resid[idxs[i]]
      end
      StatsPlots.@series begin
        seriestype := :histogram
        xaxis --> "Residuals"
        yaxis --> "Frequency"
        title --> "Histogram of Residuals"
        label --> "Observed Values"
        linewidth := 0.7
        normalize := :pdf
        ylims := (0, Inf)
        bar_position := :stack
        subplot := 3
        label := labels[i]
        color := color
        resid[idxs[i]]
      end
      @series begin
        seriestype := :scatter
        label := labels[i]
        color := color
        subplot := 4
        theoritical_quantiles, empirical_quantiles
      end
    end
  end
  @series begin
    seriestype := :hline
    xaxis --> string("Fitted ", y_name)
    yaxis --> "Residuals"
    title --> "Residuals vs Fitted"
    label --> :none
    seriesalpha := 1.0
    linecolor := :darkgrey
    line --> 2
    subplot := 2
    [0], [0]
  end
  @series begin
    seriestype := :line
    xaxis --> "Theoretical Quantiles"
    yaxis --> "Empirical Quantiles"
    title --> "Normal Q-Q"
    label --> "Q-Q Line"
    seriesalpha := 1.0
    line := 2
    linecolor := :darkgrey
    subplot := 4
    qqline, qqline
  end
end

@inline graph(model::FittedLinearModel) = plot(model::FittedLinearModel)