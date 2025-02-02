# Getting Started

## Computing the Cubage

The [`cubage`](@ref) function calculates the volume of a tree (cubage) using different cubing methods. It can handle both single trees and multiple trees in a dataset. The available cubing methods are [`Smalian`](@ref), [`Newton`](@ref) and [`Huber`](@ref).

### Cubing a Simple Tree

When cubing a single tree, you need to provide vectors of diameters (d) and heights (h) measured at different points along the tree stem. Diameters should be in centimeters, and heights should be in meters. The diameter at breast height (DBH) and total tree height (Ht) are essential inputs.

```@example ex_cub_01
using ForestMensuration

# Diameters at different heights (cm)
d = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0]

# Corresponding heights (m)
h =  [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8]

# Calculate cubage using the Smalian method
cubage(Smalian, h, d)
```

\

- vt: Total volume
- v0: Volume of the stump
- vc: Commercial bole volume
- vr: Residual volume above commercial limit
- vn: Volume of the top (cone)
- dbh: Diameter at breast height
- ht: Total tree height
- hc: Commercial height
- aff: Artificial form factor
- nff: Natural form factor
- qf: Form quotient

\

#### Including Bark Thickness

With the bark thickness value, it is possible to calculate the bark factor and total and commercial volumes without bark. Note: the provided thickness should be the 'single thickness' in centimeters. The function will convert it into 'double thickness'.

```@example ex_cub_01
# Bark thickness at corresponding heights (cm)
bark = [0.9, 0.5, 0.3, 0.2, 0.2, 0.1, 0.0]

# Define a commercial diameter limit
diameter_limit = 4.0

# Calculate cubage using the Newton method, including bark thickness and diameter limit
cubage(Newton, h, d, bark, diameter_limit)
```

\

Additional columns include:

- k: Bark factor
- : vtwb, v0wb, vcwb, vrwb, vnwb: Corresponding volumes without bark

### Cubing Multiple Trees

To calculate cubage for multiple trees, organize your data in a DataFrame with columns for tree identifiers, heights, diameters, and optionally bark thickness.

```@example ex_cub_02
using ForestMensuration
using DataFrames

# Sample data for multiple trees
data = DataFrame(
    tree = [148, 148, 148, 148, 148, 148, 148, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222],
    h = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8, 0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 11.3, 13.3, 15.3, 17.3, 19.5],
    d = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0, 16.0, 12.0, 11.6, 10.1, 9.4, 8.2, 7.0, 6.0, 4.0, 2.0, 0.0],
    bark = [0.9, 0.5, 0.3, 0.2, 0.2, 0.1, 0.0, 1.2, 0.5, 0.3, 0.3, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0, 0.0]
)

# Define a commercial diameter limit
diameter_limit = 2.5

# Calculate cubage for each tree using the Huber method
cubage(Huber, :tree, :h, :d, data, diameter_limit)
```

\

#### Including Bark Thickness for Multiple Trees

Additionally, bark thickness values can be provided to calculate bark factors and volumes without bark.

```@example ex_cub_02
# Calculate cubage including bark thickness
cubage(Huber, :tree, :h, :d, :bark, data, diameter_limit)
```

## Fitting Linear Regressions

The [`regression`](@ref) function automatically generates and evaluates multiple regression models based on the provided data. It explores various transformations of the dependent and independent variables, creating a comprehensive set of models for analysis.

### Adjusting a Hypsometric Relationship

In forestry, hypsometric relationships model the relationship between tree height (h) and diameter at breast height (dbh). The regression function generates numerous models to find the best fit.

```@example regression_data
using ForestMensuration
using DataFrames

# Sample dataset with tree heights and diameters
data = DataFrame(
    plot = ["A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",
            "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B",
            "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C",
            "D", "D", "D", "D", "D", "D", "D", "D", "D", "D", "D", "D", "D", "D", "D"],
    h = [20.9, 19.6, 13.2, 23.3, 19.2, 16.2, 8.3, 19.7, 11.0, 24.0, 25.8, 28.2, 24.2, 26.2, 28.3,
         14.4, 14.9, 15.6, 8.2, 22.1, 16.7, 22.3, 19.5, 15.9, 16.7, 24.5, 21.7, 23.8,
         20.8, 17.7, 19.3, 16.7, 22.2, 18.6, 6.9, 22.3, 8.7, 22.1, 21.0, 23.5,
         19.5, 19.7, 18.2, 13.9, 12.3, 14.5, 12.3, 18.6, 18.0, 17.4, 24.3, 22.8, 23.2, 23.5, 25.2],
    dbh = [31.5, 30.0, 26.5, 31.0, 29.0, 26.5, 14.5, 28.8, 19.0, 31.5, 32.5, 33.8, 32.5, 33.3, 36.0,
           24.0, 28.0, 23.0, 15.5, 31.0, 27.0, 29.0, 28.0, 26.0, 29.0, 30.0, 29.0, 30.5,
           25.0, 26.8, 27.5, 26.0, 26.0, 25.8, 10.8, 27.0, 16.5, 26.5, 27.0, 26.3,
           26.0, 25.5, 25.0, 23.5, 22.0, 23.0, 23.0, 26.0, 25.5, 27.5, 26.5, 26.5, 27.8, 26.0, 27.0]
)

# Perform regression analysis between height and diameter
models = regression(:h, :dbh, data)
```

This generates 240 different models combining various transformations of h and dbh.

```@example regression_data
#number of fitted regressions
length(models)
```

#### Regression Selection Criteria

After fitting the models, you can evaluate and rank them based on specific criteria using the [`criteria_table`](@ref) function.

```@example regression_data
# Evaluate models
best_models = criteria_table(models)
```

```@example regression_data
# Evaluate models based on Adjusted R², Standard Error and chosing the 5 bests
best_5_models = criteria_table(models, :adjr2, :syx, best=5)
```

\

#### Selecting the Best Model

To select the best model based on the combined ranking you can simply use the [`criteria_selection`](@ref) function:

```@example regression_data
# Select the top model
top_model = criteria_selection(models)
```

#### Plotting the Regression

You can visualize the regression model using the [`plot_regression`](@ref) function.

```@example regression_data
# Plot the top model
plot_regression(top_model)
```

\

#### Predict

The [`predict`](@ref) function. allows you to generate predicted values from a regression model on the original scale of the dependent variable. This is particularly useful when the model involves transformations of the dependent variable (e.g., logarithmic transformations). The function automatically applies the appropriate inverse transformations and corrections, such as the Meyer correction factor for logarithmic models.

```@example regression_data
# Returns the predicted values from the model on the original scale
h_pred = predict(top_model)
```

The [`predict!`](@ref) function extends this by adding the predicted values directly to your DataFrame. It creates new columns for the predicted and actual values, combining observed measurements with model predictions where data may be missing. This is especially useful in forest inventory datasets where certain tree attributes might not be measured for every tree, and predictions need to be filled in for these gaps.

```@example regression_data
# Automatically adds predicted and actual height columns to the provided DataFrame.
# This combines observed heights and predicted heights for trees with missing or unmeasured heights.
predict!(top_model, data)

# Firsts values of dataset
println(data[1:10, :])
```

### Adjusting a Qualitative (Dummy) Hypsometric Relationship

If your data includes categorical variables (e.g., different plots or species), you can include them in the regression analysis.

```@example regression_data
# Perform regression including 'plot' as a categorical variable
qualitative_models = regression(:h, :dbh, data, :plot)

# Select the best model
top_qual_model = criteria_selection(qualitative_models, :adjr2, :syx, :aic)
```

\

#### Plotting the Qualitative Regression

```@example regression_data
# Plot the top qualy model
plot_regression(top_qual_model)
```

## Site Classification

## Frequency and Statistical Functions

### Creating Frequency Tables

The [`frequency_table`](@ref) function creates frequency distributions for a vector of values, which is useful for analyzing the distribution of diameters or heights in your data.

```@example regression_data
# Frequency table for diameters using Sturges' formula for class intervals
frequency_table(data.dbh)
```

\

- LI: Lower class limit
- Xi: Class center
- LS: Upper class limit
- fi: Frequency count
- Fi: Cumulative frequency
- fri: Relative frequency (%)
- Fri: Cumulative relative frequency (%)

\

#### Specifying Class Width

You can specify the class width (hi) to customize the intervals.

```@example regression_data
# Frequency table for heights with class width of 4 meters
frequency_table(data.h, 4)
```

## Calculating Dendrometric Averages

he [`dendrometric_averages`](@ref) function computes various dendrometric averages of a forest stand, providing insights into the stand structure and growth patterns.

```@example regression_data
# Calculate dendrometric averages for the dataset
dendrometric_averages(data.dbh, area=0.05)
```

\

- d₋: Lower Hohenadl's diameter
- d̄: Mean diameter
- dg: Quadratic mean diameter
- dw: Weise's diameter (60th percentile)
- dz: Diameter of the tree with central basal area
- d₁₀₀: Mean diameter of the 100 largest trees per hectare (returns NaN if fewer than 100 trees)
- d₊: Upper Hohenadl's diameter

\

#### Estimating Heights Using a Regression Model

If you have a regression model (e.g., top_model from earlier), you can estimate the corresponding heights for each calculated diameter.

```@example regression_data
# Estimate heights for dendrometric averages using the regression model
dendrometric_averages(top_model, area=0.05)
```
