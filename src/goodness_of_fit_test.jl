"""
    goodness_of_fit_test(dist::ContinuousUnivariateDistribution, x::Vector{<:Real})

Performs a goodness-of-fit test for the given distribution `dist` against the observed data `x`. The function automatically selects the appropriate statistical test based on the characteristics of the data:

- **Anderson-Darling Test**: Used when there are tied values in `x`.
- **Kolmogorov-Smirnov Test**:
    - **Exact Test**: Applied when all values in `x` are unique and the sample size `n` is less than 40.
    - **Approximate Test**: Applied when all values in `x` are unique and the sample size `n` is 40 or greater.

After selecting and performing the appropriate test, the function calculates various evaluation metrics, including p-value, AIC, BIC, RMSE, RME, MAE, and RSE.

# Arguments
- `dist`: The fitted univariate distribution (e.g., `Normal`).
- `x`: A vector of observed data points.

# Returns
- `test`: The goodness-of-fit test.
"""
function goodness_of_fit_test(dist::ContinuousUnivariateDistribution, x::Vector{<:Real})
  n = length(x)
  # Determine the appropriate goodness-of-fit test
  if !allunique(x)
    # Use Anderson-Darling test when there are tied values
    test = OneSampleADTest(x, dist)
  else
    if n < 100
      # Use Exact Kolmogorov-Smirnov test for small sample sizes
      test = ExactOneSampleKSTest(x, dist)
    else
      # Use Approximate Kolmogorov-Smirnov test for larger sample sizes
      test = ApproximateOneSampleKSTest(x, dist)
    end
  end

  return test
end