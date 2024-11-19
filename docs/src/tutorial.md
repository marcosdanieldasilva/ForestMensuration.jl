# Getting Started

## Computing the [`cubage`](@ref)

### Cubing a Simple Tree

Cubage for a single tree, with diameters measured in centimeters and heights in meters. The diameter at breast height (DBH) and total height (Ht) must be provided as part of the input values.
The CubingMethods can be: [`Smalian`](@ref), [`Newton`](@ref) and [`Huber`](@ref).

```@example ex_cub_01
using ForestMensuration

d = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0]

h =  [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8]

cubage(Smalian, h, d)
```

With the bark thickness value, it is possible to calculate the bark factor and total and commercial volumes without bark. Note: the provided thickness should be the 'single thickness' in millimeters, i.e., the actual value collected from the field. The function will convert it into 'double thickness'.

```@example ex_cub_01

bark = [0.9, 0.5, 0.3, 0.2, 0.2, 0.1, 0.0]

cubage(Newton, h, d, bark)
```

### Cubing Multiple Trees

We utilize the function cubage to calculate cubage for multiple trees. The function expects data in a dataframe format where each row represents a tree, and columns contain attributes such as height and diameter. The cubage function can be applied to each tree group using the specified height and diameter columns.

```@example ex_cub_02
using ForestMensuration # hide
using DataFrames

data = DataFrame(
    tree = [148, 148, 148, 148, 148, 148, 148, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222],
    h = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8, 0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 11.3, 13.3, 15.3, 17.3, 19.5],
    d = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0, 16.0, 12.0, 11.6, 10.1, 9.4, 8.2, 7.0, 6.0, 4.0, 2.0, 0.0],
    bark = [0.9, 0.5, 0.3, 0.2, 0.2, 0.1, 0.0, 1.2, 0.5, 0.3, 0.3, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0, 0.0]
)

cubage(Huber, :tree, :h, :d, data)
```

Additionally, bark thickness values can be provided to calculate bark factors and volumes without bark.

```@example ex_cub_02
cubage(Huber, :tree, :h, :d, :bark, data)
```

## Fitting the [`regression`](@ref)

### Adjusting a Hypsometric Relationship

The regression function in ForestMensuration.jl is designed to explore multiple transformations and model combinations between two quantitative variables. In forestry, this function can be particularly useful for assessing relationships such as height (h) in function of diameter at breast height (dbh), as it automatically generates 240 unique combinations of transformations for both dependent and independent variables.

By analyzing a wide range of transformations, including logarithmic, inverse, and square root forms, the regression function allows for a comprehensive assessment of potential relationships in your data, aiding in the selection of the best-fit model.

```@example regression_data
using ForestMensuration
using DataFrames

# Full dataset with 100 values, 3 plots
data = DataFrame(
    h = [35.2, 34.5, 28.8, 28.5, 33.5, 29.5, 28.4, 30.4, 33.0, 31.0,
         24.8, 36.9, 34.0, 28.2, 28.9, 31.5, 28.0, 29.5, 34.6, 29.6,
    ],
    dbh = [42.0, 25.0, 33.5, 21.3, 27.0, 24.0, 23.3, 39.0, 42.0, 42.5,
           23.0, 46.0, 39.0, 34.0, 26.5, 41.0, 31.3, 30.5, 30.5, 28.5
    ]
)

# Performing the regression analysis on the full dataset
# Here, we analyze the relationship between height (h) and diameter (dbh)
models = regression(:h, :dbh, data);
# Alternative print of fitted models
models_eq = ModelEquation.(models)
```

```@example regression_data
#number of fitted regressions
length(models)
```

### Viewing the top models based on specific criteria

```@example regression_data
# Using all criteria and presenting the 10 best models
best_models = criteria_table(models)
```

```@example regression_data
# Chossing as especific criteria and best 5, if best = false or 0 show all regressions from 'models'
best_models_v2 = criteria_table(models, :adjr2, :rmse, best=5)
```

### Regression Selection Criteria [`criteria_table`](@ref)

### Site Class

## Frequency and Diametric Tables

## Inventory Report
