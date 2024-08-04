# function _confidence_interval(fitted_model::FittedLinearModel)
#   model_matrix = modelmatrix(fitted_model)
#   R = fitted_model.chol.U
#   variance = deviance(fitted_model) / dof_residual(fitted_model)
#   residvar = ones(size(model_matrix, 2)) * variance
#   tvalue = -quantile(TDist(dof_residual(fitted_model)), 0.025)
#   confidence = tvalue * .√((model_matrix / R).^2 * residvar)
#   return  confidence
# end

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
      color = cgrad(:darktest, ngroups, categorical = true)[i]
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

@inline graph(fitted_model::FittedLinearModel; kw...) = plot(fitted_model; kw...)

@recipe function f(fitted_model::FittedLinearModel, split::Bool)
  if !split
      pred = predict(fitted_model)
      resid = residuals(fitted_model)
      confidence = confidence_interval(fitted_model)
      y_name, x_name = names(fitted_model.data)[1:2]
      n = size(fitted_model.data, 1)
      qp = qqbuild(fit_mle(Normal, resid), resid)
      x_qqline = [extrema(qp.qx)...]
      quantx, quanty = quantile(qp.qx, [0.25, 0.75]), quantile(qp.qy, [0.25, 0.75])
      slp = diff(quanty) ./ diff(quantx)
      y_qqline = quanty .+ slp .* (x_qqline .- quantx)
      grid --> :none
      tick_direction --> :out
      framestyle --> :box
      markerstrokewidth --> 0
      seriesalpha --> 0.6
      legend --> :outertopright
      (y, x, q...) = eachcol(fitted_model.data)
      color_palette --> cgrad(:darktest, categorical = true)
      # titlefontize --> 10
      guidefontsize --> 9
      legendfontsize --> 7
      fontfamily --> :times
      size --> (1000, 650)
      margin --> 3mm
      layout := (2, 2)
      if size(fitted_model.data, 2) == 2
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
              seriestype := :hline
              subplot := 2
              label := :none
              seriesalpha := 1.0
              linecolor := :darkgrey
              line --> 2
              [0], [0]
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
              qp.qx, qp.qy
          end
          @series begin
              seriestype := :line
              subplot := 4
              label --> :none
              seriesalpha := 1.0
              line := 2
              linecolor := :darkgrey
              x_qqline, y_qqline
          end
      else
          h = Plots._make_hist((vec(copy(resid)),), :auto; normed = true)
          nbins = length(h.weights)
          edges = h.edges[1]
          bar_width = mean(map(i -> edges[i + 1] - edges[i], 1:nbins))
          (y, x, q...) = eachcol(fitted_model.data)
          x_bar = map(i -> (edges[i] + edges[i + 1]) / 2, 1:nbins)
          gb = RecipesPipeline._extract_group_attributes(tuple(q...))
          labels, idxs = getfield(gb, 1), getfield(gb, 2)
          ngroups = length(labels)
          ntot = count(x_bar -> !isnan(x_bar), resid)
          y_bar = fill(NaN, nbins, ngroups)
          for i = 1:ngroups
              groupinds = idxs[i]
              v_i = filter(x_bar -> !isnan(x_bar), resid[groupinds])
              h_i = Plots._make_hist((v_i,), h.edges; normed = false, weights = nothing)
              y_bar[:, i] .= h_i.weights .* (length(v_i) / ntot / sum(h_i.weights))
              y_bar, fr = StatsPlots.groupedbar_fillrange(y_bar)
          end
          for i = 1:ngroups
              colors = cgrad(:darktest, categorical = true, ngroups)
              qp = qqbuild(fit_mle(Normal, resid[idxs[i]]), resid[idxs[i]])
              @series begin
                  seriestype := :scatter
                  subplot := 1
                  label := labels[i]
                  color := colors[i]
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
                  color := colors[i]
                  xvals = x[idxs[i]]
                  predvals = pred[idxs[i]]
                  xvals, predvals
              end
              @series begin
                  seriestype := :scatter
                  subplot := 2
                  label := labels[i]
                  color := colors[i]
                  pred[idxs[i]], resid[idxs[i]]
              end
              @series begin
                  seriestype := :bar
                  seriesalpha := 1.0
                  subplot := 3
                  label := labels[i]
                  color := colors[i]
                  yvals = y_bar[:, i]
                  x_bar, yvals
              end
              @series begin
                  seriestype := :scatter
                  label := labels[i]
                  color := colors[i]
                  subplot := 4
                  qp.qx, qp.qy
              end
          end
      end
  else
      pred = predict(fitted_model)
      resid = residuals(fitted_model)
      confidence = confidence_interval(fitted_model)
      (y, x, q...) = eachcol(fitted_model.data)
      y_name, x_name = names(fitted_model.data)[1:2]
      gb = RecipesPipeline._extract_group_attributes(tuple(q...))
      labels, idxs = getfield(gb, 1), getfield(gb, 2)
      ngroups = length(labels)
      colors = cgrad(:darktest, categorical = true, ngroups)
      layout := (4, ngroups)
      grid --> :none
      tick_direction --> :out
      framestyle --> :box
      markerstrokewidth --> 0
      seriesalpha --> 0.6
      for i = 1:ngroups
          qp = qqbuild(fit_mle(Normal, resid[idxs[i]]), resid[idxs[i]])
          x_qqline = [extrema(qp.qx)...]
          quantx, quanty = quantile(qp.qx, [0.25, 0.75]), quantile(qp.qy, [0.25, 0.75])
          slp = diff(quanty) ./ diff(quantx)
          y_qqline = quanty .+ slp .* (x_qqline .- quantx)
          @series begin
              seriestype := :scatter
              subplot := i
              link := :y
              label := :none
              color := colors[i]
              xvals = x[idxs[i]]
              yvals = y[idxs[i]]
              xvals, yvals
          end
          @series begin
              seriestype := :line
              subplot := i
              ribbon := confidence[idxs[i]]
              seriesalpha := 0.7
              line --> 2
              title --> labels[i]
              label := :none
              color := colors[i]
              xvals = x[idxs[i]]
              predvals = pred[idxs[i]]
              xvals, predvals
          end
          @series begin
              seriestype := :scatter
              subplot := i + ngroups
              label := :none
              # link := :both
              color := colors[i]
              pred[idxs[i]], resid[idxs[i]]
          end
          @series begin
              seriestype := :hline
              subplot := i + ngroups
              label := :none
              seriesalpha := 1.0
              linecolor := :darkgrey
              line --> 2
              [0], [0]
          end
          @series begin
              seriestype := :histogram
              subplot := i + 2ngroups
              label := :none
              linewidth := 0.7
              normalize := :pdf
              color := colors[i]
              ylims := (0, Inf)
              resid[idxs[i]]
          end
          @series begin
              seriestype := :scatter
              subplot := i + 3ngroups
              label := :none
              color := colors[i]
              qp.qx, qp.qy
          end
          @series begin
              seriestype := :line
              subplot := i + 3ngroups
              label --> :none
              seriesalpha := 1.0
              line := 2
              linecolor := :darkgrey
              x_qqline, y_qqline
          end
      end
  end
end

@inline graph(fitted_model::FittedLinearModel, split::Bool=false; kw...) = Plots.plot(fitted_model, split; kw...)
