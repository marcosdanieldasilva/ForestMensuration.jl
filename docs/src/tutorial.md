# Getting Started

## Computing the Cubage

The `cubage` function calculates the volume of a tree by sections from a profile of diameters and heights. It partitions the volume into stump, commercial, residual, and top sections, and can handle both single and multiple trees. The volume of each section is calculated using Smalian's formula. For calculating the volume of individual logs with other methods, see `smalian`, `newton`, and `huber`.

### Cubing a Simple Tree

When cubing a single tree, you need to provide vectors of diameters (d) and heights (h) measured at different points along the tree stem. Diameters should be in centimeters, and heights should be in meters. The diameter at breast height (DBH) and total tree height (Ht) are essential inputs.

```@example ex_cub_01
using ForestMensuration
using Unitful

# Diameters at different heights (cm)
d = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0]u"cm"

# Corresponding heights (m)
h =  [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8]u"m"

# Calculate cubage
cubage(h, d)
```

\

- vt: Total volume
- v0: Volume of the stump
- vc: Commercial bole volume
- vr: Residual volume above commercial limit
- vn: Volume of the top (cone)
- d: Diameter at breast height
- h: Total tree height
- hc: Commercial height
- aff: Artificial form factor
- nff: Natural form factor
- qf: Form quotient

\

#### Including Bark Thickness

With the bark thickness value, it is possible to calculate the bark factor and total and commercial volumes without bark. Note: the provided thickness should be the 'single thickness' in centimeters. The function will convert it into 'double thickness'.

```@example ex_cub_01
# Bark thickness at corresponding heights (cm)
bark = [0.9, 0.5, 0.3, 0.2, 0.2, 0.1, 0.0]u"cm"

# Define a commercial diameter limit
diameter_limit = 4.0u"cm"

# Calculate cubage including bark thickness and diameter limit
cubage(h, d, bark; dlimit = diameter_limit)
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
using Unitful

# Sample data for multiple trees
data = DataFrame(
    tree = [148, 148, 148, 148, 148, 148, 148, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222],
    h = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8, 0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 11.3, 13.3, 15.3, 17.3, 19.5],
    d = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0, 16.0, 12.0, 11.6, 10.1, 9.4, 8.2, 7.0, 6.0, 4.0, 2.0, 0.0],
    bark = [0.9, 0.5, 0.3, 0.2, 0.2, 0.1, 0.0, 1.2, 0.5, 0.3, 0.3, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0, 0.0]
)

# Define a commercial diameter limit
diameter_limit = 2.5u"cm"

# Calculate cubage for each tree
cubage(data.tree, data.h .* u"m", data.d .* u"cm"; dlimit = diameter_limit)
```

\

#### Including Bark Thickness for Multiple Trees

Additionally, bark thickness values can be provided to calculate bark factors and volumes without bark.

```@example ex_cub_02
# Calculate cubage including bark thickness
cubage(data.tree, data.h .* u"m", data.d .* u"cm", data.bark .* u"cm"; dlimit = diameter_limit)
```

## Fitting Allometric Regressions

The [`regression`](@ref) function — from
[ForestModeling.jl](https://github.com/JuliaForests/ForestModeling.jl), re-exported here —
automatically generates and evaluates multiple regression models based on the provided
data. It explores a bounded, citable catalog of transformations of the dependent and
independent variables, creating a comprehensive set of candidate models for analysis.

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

# Perform regression analysis between height and diameter (data comes first, then y, then x...)
models = regression(data, :h, :dbh)
```

This generates every combination (up to `nMax=2` terms, the default) of the transform
catalog applied to h and dbh.

```@example regression_data
#number of fitted regressions
length(models)
```

#### Regression Selection Criteria

After fitting the models, you can evaluate and rank them based on specific criteria using the [`criteriaTable`](@ref) function.

```@example regression_data
# Evaluate models
best_models = criteriaTable(models)
```

```@example regression_data
# Evaluate models based on Adjusted R², Coefficient of Variation, choosing the 5 best
best_5_models = criteriaTable(models, :adjr2, :cv; best=5)
```

\

#### Selecting the Best Model

To select the best model based on the combined ranking you can simply use the [`criteriaSelection`](@ref) function:

```@example regression_data
# Select the top model
top_model = criteriaSelection(models, :adjr2, :cv)
```

#### Predict

The [`predict`](@ref) function allows you to generate predicted values from a regression model on the original scale of the dependent variable. This is particularly useful when the model involves transformations of the dependent variable (e.g., logarithmic transformations). The function automatically applies the appropriate inverse transformations and bias corrections.

```@example regression_data
# Returns the predicted values from the model on the original scale
h_pred = predict(top_model)
```

### Adjusting a Qualitative (Dummy) Hypsometric Relationship

If your data includes categorical variables (e.g., different plots or species), you can include them in the regression analysis simply by passing them alongside the continuous predictors — the engine detects categorical columns from the data automatically.

```@example regression_data
# Perform regression including 'plot' as a categorical variable
qualitative_models = regression(data, :h, :dbh, :plot)

# Select the best model
top_qual_model = criteriaSelection(qualitative_models, :adjr2, :cv)
```

### Robust Regression (Alternative Estimation Criteria)

The [`fitRobust`](@ref) function fits a single formula by minimizing an alternative loss
with `Optim.jl` instead of closed-form least squares — useful when OLS assumptions are
visibly violated (heavy outliers, or strongly relative/percentage-scale error). It accepts
`SSE` (reproduces OLS), `MAE`, `HUBER`, `MSLE`, or `MAPE`.

```@example regression_data
# Huber loss down-weights outliers relative to ordinary least squares
robust_model = fitRobust(@formula(log(h) ~ log(dbh)), data, HUBER)
```

```@example regression_data
predict(robust_model)
```

`fitRobust` is opt-in — [`regression`](@ref) never calls it, since numerical optimization
is far more expensive to run combinatorially than the closed-form path.

### Regression with Units of Measurement

Every regression function accepts `Unitful` quantity columns directly, exactly like the
rest of ForestMensuration.jl — units are stripped before fitting and reattached to
`predict`/`fitted`/`residuals`, since the transform catalog (`log`, `1/x`, `√x`, ...)
cannot run on a dimensioned quantity.

```@example regression_units
using ForestMensuration, DataFrames, Unitful

data_u = DataFrame(
    dbh = [31.5, 30.0, 26.5, 31.0, 29.0, 26.5, 14.5, 28.8, 19.0, 31.5]u"cm",
    h   = [20.9, 19.6, 13.2, 23.3, 19.2, 16.2, 8.3, 19.7, 11.0, 24.0]u"m",
)

models_u = regression(data_u, :h, :dbh; nMax=2)
best_u = criteriaSelection(models_u, :adjr2, :cv)
```

```@example regression_units
# predictions come back as Unitful quantities, on the same scale h was fit in
predict(best_u)
```

```@example regression_units
# a model fit on cm/m scores a table given in a compatible unit automatically
predict(best_u, DataFrame(dbh=[300.0]u"mm"))
```

Summary statistics and tables (`rmse`, `mae`, `criteriaTable`, `metrics`, `siteTable`) stay
plain numbers regardless — a coefficient built from an arbitrarily transformed predictor
has no single clean physical unit.

## Grouped/Stratified Regression

For stratified data (e.g. several species that plausibly need different equations),
ForestModeling.jl — re-exported here — provides a grouped/stratified variant of the same
bounded transform search: the [`regressionGrouped`](@ref) function fits three strategies
so you can tell whether stratifying is actually worth it — one pooled equation
(`general`), one pooled equation with the group as a categorical covariate (`qualy`), and
one independently selected equation per group (`grouped`).

```@example regression_grouped
using ForestMensuration, DataFrames

# 3 species with distinct slope/intercept, so stratifying is expected to help
d = repeat(10.0:3.0:37.0, 3)
species = repeat(["Oak", "Pine", "Cedar"], inner=10)
slope = Dict("Oak" => 0.85, "Pine" => 0.65, "Cedar" => 0.75)
intercept = Dict("Oak" => 0.2, "Pine" => 0.5, "Cedar" => 0.35)
h = [intercept[species[i]] + slope[species[i]] * log(d[i]) + 0.03 * sin(i) for i in 1:30]
gdata = DataFrame(d=d, h=h, species=species)

grouped_model = regressionGrouped(gdata, :h, :d, :species)
criteriaTable(grouped_model, :adjr2, :cv)
```

```@example regression_grouped
# per-group breakdown instead of the pooled general/qualy/grouped comparison
criteriaTable(collect(values(grouped_model.grouped)), :adjr2)
```

### Range-Safe Prediction

Calling [`predict`](@ref) directly on a small group's own equation risks wild
extrapolation the moment a new row's predictor falls outside that group's own fitted
range. [`predictBounded`](@ref) checks the range per row instead: out of the **global**
range → the pooled `general` model; in the global range but out of that row's **own
group's** range (or no model for that group at all) → `qualy`; only otherwise does the
group's own model get used — so a small subgroup can't extrapolate into nonsensical
predictions.

```@example regression_grouped
new_trees = DataFrame(d=[15.0, 60.0, 15.0], species=["Cedar", "Cedar", "Birch"])
predictBounded(grouped_model, new_trees)
```

## Site Classification

Site classification lets you evaluate and classify forest sites based on regression
models relating tree height and age. These functions are particularly useful for
assessing site productivity and quality by comparing observed data with expected values
derived from a well-calibrated model.

### Calculating Site Classification

The [`siteClassification`](@ref) function calculates the expected dominant height at a given index age for each observation based on a fitted regression model. This is a key step in classifying the productivity of a forest site. As with every function in this package, heights are given as `Unitful` quantities (`u"m"` below) so the result carries the same unit automatically — plain numbers work too and are assumed to already be in meters, but the unitful version is the recommended way to call it.

```@example site_classification
using ForestMensuration, DataFrames

# Create a DataFrame containing tree plot data -- height in meters
data = DataFrame(
    plot = repeat(1:6, inner=5),
    age  = repeat([36, 48, 60, 72, 84], outer=6),
    h    = [13.6, 17.8, 21.5, 21.5, 21.8,
            14.3, 17.8, 21.0, 21.0, 21.4,
            14.0, 17.5, 21.2, 21.2, 21.4,
            13.4, 18.0, 20.8, 20.8, 23.2,
            13.2, 17.4, 20.3, 20.3, 22.0,
            13.2, 17.8, 21.3, 21.3, 22.5]u"m"
)

# Fit a regression model to relate height (h) to age, and pick the best one
reg = criteriaSelection(regression(data, :h, :age), :adjr2, :cv)

# Define the target index age (for example, 60 months)
index_age = 60

# Calculate the site classification values (site indices) for each observation -- a Vector of Unitful heights
site_indices = siteClassification(reg, data, index_age)

println("Site Classification Values:")
println(site_indices)
```

### Calculating Dominant Height Classification

The [`hdomClassification`](@ref) function uses the site classification values to predict the dominant height for each observation at the specified index age. This reverses the site classification process, allowing you to forecast tree heights based on site productivity.

```@example site_classification
# Now, compute the dominant heights for each observation using the site indices
dominant_heights = hdomClassification(reg, data, index_age, site_indices)

println("Dominant Height Values:")
println(dominant_heights)
```

### Generating a Site Table

The [`siteTable`](@ref) function creates a table of predicted dominant heights at various
ages for different site index classes. You can specify a height increment (`hi`) to
define the granularity of the site classes; it is chosen automatically via Sturges' rule
when omitted. Unlike `siteClassification`/`hdomClassification` above, this table's values
stay plain `Float64` in the fitting unit (meters here) rather than `Unitful` quantities —
the same convention `criteriaTable`/`metrics` use in ForestModeling.jl for summary tables.

```@example site_classification
# Generate the site table -- values are plain numbers (meters), by design
site_table = siteTable(reg, index_age)
```

## Frequency and Statistical Functions

### Creating Frequency Tables

The [`frequencytable`](@ref) function creates frequency distributions for a vector of values, which is useful for analyzing the distribution of diameters or heights in your data.

```@example regression_data
# Frequency table for diameters using Sturges' formula for class intervals
frequencytable(data.dbh)
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
frequencytable(data.h, 4)
```

## Calculating Dendrometric Averages

To evaluate the horizontal and vertical structure of a forest stand, ForestMensuration.jl provides comprehensive metric summary functions. We can analyze diameters, heights, and complete plot-level inventories in a single step.

First, let's define our sample plot data using Unitful quantities:

```@example metrics
using ForestMensuration

diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0] * u"cm"
heights = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2] * u"m"
plotArea = 500u"m^2"
nothing # hide

```

### Horizontal Structure (Diameters)

The [`dmetrics`](@ref) function calculates a comprehensive set of dendrometric diameter averages, returning a single-row DataFrame.

```@example metrics
# Calculate diameter metrics
dmetrics(diameters, plotArea)

```

The resulting columns represent:

- **dl**: Lower Hohenadl diameter (mean minus standard deviation).
- **dm**: Arithmetic mean diameter.
- **dg**: Quadratic mean diameter (diameter of the tree with mean basal area).
- **dw**: Weise's diameter (60th percentile).
- **dz**: Central basal area diameter.
- **dd**: Dominant diameter (mean of the thickest trees).
- **du**: Upper Hohenadl diameter (mean plus standard deviation).
- **dv**: Coefficient of variation of the diameters (%).

### Vertical Structure (Heights)

Similarly, the [`hmetrics`](@ref) function calculates the vertical structure metrics of the stand. It requires both diameters and heights to calculate weighted averages (like Lorey's mean height) and Assmann's dominant height.

```@example metrics
# Calculate height metrics
hmetrics(diameters, heights, plotArea)

```

The resulting columns represent:

- **hl**: Lower height boundary (mean minus standard deviation).
- **hm**: Arithmetic mean height.
- **hd**: Dominant height (based on the 100 thickest trees per hectare).
- **hg**: Lorey's mean height (weighted by basal area).
- **hu**: Upper height boundary (mean plus standard deviation).
- **hv**: Coefficient of variation of the heights (%).

### Comprehensive Stand Summary

For a complete characterization of an individual plot, use the [`standmetrics`](@ref) function. It unifies horizontal structure, vertical structure, and spatial density into a single output, automatically scaling plot totals to per-hectare or per-acre equivalents using an Expansion Factor (EF).

```@example metrics
# Calculate full stand metrics and per-area extrapolations
standmetrics(diameters, heights, plotArea)

```

In addition to all the structural variables from `dmetrics` and `hmetrics`, the summary adds:

- **n**: Total number of sampled trees in the plot.
- **g**: Total basal area of the sampled plot.
- **EF**: Expansion factor used to scale plot data to a per-hectare (or per-acre) basis.
- **N**: Extrapolated number of trees per hectare/acre.
- **G**: Extrapolated total basal area per hectare/acre.

## Forest Inventory Sampling

ForestMensuration.jl implements all 11 classic forest inventory sampling designs, each
generic to any number of strata/clusters/plots. Every design accepts plain numbers
(volume defaults to `m^3`, plot/total areas to `ha`) or explicit `Unitful` quantities.
Designs that produce a single table return a plain `DataFrame`; designs that produce
several related tables (one per stratum/cluster, or one per inventory occasion) return a
[`SamplingReport`](@ref), whose tables are reachable both as `report.tables.name` and
directly as `report.name`.

### Simple Random Sampling

The [`simplecasualsampling`](@ref) function is the reference design every other method in
this section is compared against: each plot has an equal chance of being selected, with
no further structure.

```@example inv_simple
using ForestMensuration

v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]

simplecasualsampling(v, 0.05, 10; e=10, α=0.95)
```

### Stratified Random Sampling

The [`stratifiedsampling`](@ref) function divides the population into non-overlapping
strata before sampling within each — usually a variance reduction over simple random
sampling of the same total size. It returns a [`SamplingReport`](@ref) with an ANOVA
table (is there really a difference between strata?), an auxiliary table (per-stratum
descriptive statistics and allocation weights), and the final result table.

```@example inv_stratified
using ForestMensuration, DataFrames

data = DataFrame(
    stratum=[1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3],
    volume=[18.2, 21.4, 19.8, 20.1, 32.5, 35.1, 30.8, 12.4, 11.9, 13.6, 12.8, 13.1],
)

# strata_area is given in the same order as the sorted strata: 1, 2, 3
report = stratifiedsampling(:stratum, :volume, 0.1, [12.0, 8.0, 20.0], data)
resultTable(report)
```

```@example inv_stratified
auxiliaryTable(report)
```

```@example inv_stratified
anova(report)
```

With a single stratum, `stratifiedsampling` reduces exactly to `simplecasualsampling` —
useful as a sanity check when strata are added incrementally to a growing dataset.

### Systematic Sampling

The [`systematicsampling`](@ref) function estimates the variance of the mean from the
method of successive differences between consecutive plots, since plots laid out at a
fixed interval tend to be more alike than a true random sample.

```@example inv_systematic
using ForestMensuration

v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]

systematicsampling(v, 0.05, 10)
```

An optional `line` vector groups plots into several independent transects, so the
difference between the last plot of one line and the first plot of the next — which
aren't actually adjacent on the ground — is excluded from the variance estimate:

```@example inv_systematic
line = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
systematicsampling(v, 0.05, 10; line=line)
```

### Systematic Sampling with Multiple Random Starts

The [`multistartsystematicsampling`](@ref) function lays out several independent
systematic lines, each starting at its own random point, turning the design into a
genuine probability sample: statistically identical to [`clustersampling`](@ref)'s
between/within decomposition, with each line playing the role of one cluster.

```@example inv_multistart
using ForestMensuration, DataFrames

data = DataFrame(
    start=repeat(1:4, inner=5),
    volume=[18.2, 19.1, 17.8, 18.9, 19.4, 22.4, 23.1, 21.9, 22.8, 23.4,
            15.1, 16.0, 14.8, 15.6, 15.9, 20.5, 21.2, 19.8, 20.6, 21.0],
)

report = multistartsystematicsampling(:start, :volume, 0.02, 15, data)
resultTable(report)
```

### One-Stage Cluster Sampling

The [`clustersampling`](@ref) function samples clusters of `M` neighboring plots
("conglomerates") instead of individual plots — cheaper to lay out in the field, at the
cost of within-cluster homogeneity typically inflating the variance of the mean.

```@example inv_cluster
using ForestMensuration, DataFrames

data = DataFrame(
    cluster=repeat(1:6, inner=4),
    volume=[18.2, 19.1, 17.8, 18.9, 22.4, 23.1, 21.9, 22.8, 15.1, 16.0, 14.8, 15.6,
            27.3, 28.1, 26.9, 27.8, 19.8, 20.5, 19.1, 20.0, 24.5, 25.2, 23.9, 24.8],
)

report = clustersampling(:cluster, :volume, 0.02, 15, data)
resultTable(report)
```

```@example inv_cluster
clusterTable(report)
```

### Horizontal Point Sampling (Bitterlich)

The [`horizontalpointsampling`](@ref) function estimates volume and basal area per
hectare from angle-count ("Bitterlich") data: at each point, every tree that is "in" for
the chosen basal area factor (`baf`) represents exactly `baf` m²/ha of basal area,
regardless of its own size — no fixed plot radius is needed at all. Point 6 below was
visited but had no "in" trees; it is not a row of `data`, but is still counted through
`npoints=6` so the zero observation is not silently dropped from the mean.

```@example inv_hps
using ForestMensuration, DataFrames

data = DataFrame(
    point=[1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5],
    diameter=[25.0, 30.0, 20.0, 28.0, 22.0, 35.0, 30.0, 25.0, 20.0, 22.0, 30.0, 28.0, 26.0],
    volume=[0.35, 0.55, 0.22, 0.48, 0.28, 0.85, 0.55, 0.35, 0.22, 0.28, 0.55, 0.48, 0.40],
)

# baf=2 m²/ha, 6 points visited (5 counted trees + 1 empty), 0.1 ha per point, 10 ha stand
report = horizontalpointsampling(:point, :diameter, :volume, 2.0, 6, 0.1, 10, data)
resultTable(report)
```

```@example inv_hps
pointTable(report)
```

### Two-Stage Sampling

The [`twostagesampling`](@ref) function draws `n` primary units (e.g. stands) from a
population of `N`, then sub-samples `m` secondary units (plots) from within each drawn
primary out of `M` possible — unlike cluster sampling, not every secondary unit inside a
drawn primary needs to be measured.

```@example inv_twostage
using ForestMensuration, DataFrames

data = DataFrame(
    primary=repeat(1:5, inner=3),
    volume=[18.2, 19.1, 17.8, 22.4, 23.1, 21.9, 15.1, 16.0, 14.8, 27.3, 28.1, 26.9, 19.8, 20.5, 19.1],
)

# N=40 possible primary units, M=6 possible secondary units per primary
report = twostagesampling(:primary, :volume, 0.02, 40, 6, data)
resultTable(report)
```

### Successive Occasions

The remaining four designs estimate growth between two inventory occasions, differing in
how much of the plot network is remeasured at the second occasion. All four return a
[`SamplingReport`](@ref) with `occasion1`, `occasion2`, and `change` tables.

#### Independent Samples

The [`independentoccasionssampling`](@ref) function samples each occasion completely
independently — the simplest design, but the least efficient at detecting growth, since
it exploits none of the natural plot-to-plot correlation between occasions.

```@example inv_independent
using ForestMensuration

v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0]
v2 = [24.1, 23.8, 22.9, 25.6, 24.0]

report = independentoccasionssampling(v1, v2, 0.05, 200, 200)
change(report)
```

#### Complete Replacement

The [`completereplacementsampling`](@ref) function remeasures the exact same plots at
both occasions, exploiting their positive correlation to tighten the growth estimate.

```@example inv_complete
using ForestMensuration

v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 24.6, 17.3]
v2 = [22.1, 25.8, 23.4, 21.0, 26.6, 20.9, 26.0, 22.5]

report = completereplacementsampling(v1, v2, 0.05, 200, 200)
change(report)
```

#### Partial Replacement

The [`partialreplacementsampling`](@ref) function is the middle ground: a matched subset
of plots is remeasured, some temporary plots are dropped, and new temporary plots are
added at the second occasion. `missing` marks a plot that wasn't measured on a given
occasion.

```@example inv_partial
using ForestMensuration

volume1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, missing, missing]
volume2 = [missing, missing, 23.4, 24.0, 26.6, 22.9, 27.5, 25.2, 21.8]

report = partialreplacementsampling(volume1, volume2, 0.05, 200)
occasion2(report)
```

```@example inv_partial
change(report)
```

#### Double Sampling

The [`doublesampling`](@ref) function measures a large first-occasion sample but only
remeasures a smaller "permanent" subset at the second occasion, estimating the rest via a
regression of the permanent subset's second-occasion volume on its first-occasion volume
— solved through the normal equations `(X'X)β = X'y`, the same transposed-design-matrix
approach used for the ANOVA in `stratifiedsampling`.

```@example inv_double
using ForestMensuration

volume1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, 17.6, 20.8, 21.9]
volume2 = [22.1, 25.8, 23.4, missing, 26.6, missing, 27.5, missing, missing, 26.0]

report = doublesampling(volume1, volume2, 0.05, 200)
occasion2(report)
```

```@example inv_double
change(report)
```

## Stem Taper Equations

ForestMensuration.jl can fit a stem taper (profile) curve to measured `(height, diameter)`
pairs along one or more trees' stems, then use it to evaluate the diameter at any height,
invert it to find the height at any diameter, integrate a volume between two heights, and
simulate cutting the stem into logs across any number of product classes. Fitting a
[`TaperFit`](@ref) is done with `fit` from
[ForestModeling.jl](https://github.com/JuliaForests/ForestModeling.jl) (re-exported here);
everything downstream — [`taperdiameter`](@ref), [`taperheight`](@ref),
[`taperedvolume`](@ref), [`logassortment`](@ref) — is native to this package.

### Fitting a Taper Model

10 classic published taper forms are available as `TaperModel` subtypes — `Kozak1969`,
`Schoepfer1966`, `Matte1949`, `Demaerschalk1972`, `Clutter1980`, `MaxBurkhart1976`,
`Johnson1911`, `Kozak1988`, `Kozak2004`, `Bi2000`. `fit` takes paired stem-scaling data —
each tree's `dbh`/total `height` repeated once per measured section, alongside the section
heights `hi` and diameters `di`:

```@example taper_fit
using ForestMensuration

dbh    = [20.0, 20.0, 20.0, 30.0, 30.0, 30.0, 25.0, 25.0, 25.0, 35.0, 35.0, 35.0]
height = [18.0, 18.0, 18.0, 22.0, 22.0, 22.0, 20.0, 20.0, 20.0, 24.0, 24.0, 24.0]
hi     = [0.3, 6.0, 14.0, 0.3, 8.0, 18.0, 0.3, 7.0, 16.0, 0.3, 9.0, 20.0]
di     = [22.4, 15.8, 6.1, 33.6, 24.2, 8.7, 27.9, 18.6, 7.4, 38.9, 27.1, 9.9]

ft = fit(Kozak1969(), dbh, height, hi, di)
r2(ft), adjr2(ft), dispersion(ft)
```

Fitting several models on the same data and ranking them with
[`criteriaTable`](@ref)/[`criteriaSelection`](@ref) — from `ForestModeling.jl`, works
unchanged on a `Vector{TaperFit}` — picks the best-fitting form:

```@example taper_fit
fits = [fit(m, dbh, height, hi, di) for m in (Kozak1969(), Schoepfer1966(), Demaerschalk1972(), Bi2000())]
criteriaTable(fits, :adjr2, :rmse)
```

### Diameter and Height Along the Stem

[`taperdiameter`](@ref) evaluates the fitted curve; [`taperheight`](@ref) inverts it —
useful for finding, say, the height at a minimum merchantable top diameter:

```@example taper_fit
taperdiameter(ft, 25.0u"cm", 20.0u"m", 7.0u"m")
```

```@example taper_fit
taperheight(ft, 25.0u"cm", 20.0u"m", 15.0u"cm")
```

Plain numbers work the same way, assuming `cm` for diameters and `m` for heights:

```@example taper_fit
taperdiameter(ft, 25.0, 20.0, [2.0, 7.0, 12.0])
```

### Volume by Integration

[`taperedvolume`](@ref) integrates the fitted cross-sectional-area profile between two
heights — the curve-fitted counterpart of [`cubage`](@ref), which instead integrates a
*measured* section profile. Omitting `hmin`/`hmax` integrates the whole stem:

```@example taper_fit
taperedvolume(ft, 25.0u"cm", 20.0u"m", 0.3u"m", 12.0u"m")
```

```@example taper_fit
taperedvolume(ft, 25.0, 20.0)
```

### Log Assortment (Sortimentos)

[`logassortment`](@ref) simulates cutting a tree's stem into logs across any number of
product classes — a generalized, bug-fixed equivalent of the R package `timbeR`'s
per-model bucking functions. `products` is a table in **cutting priority order** (most
valuable first), one row per product, with its minimum small-end diameter (`sed`),
usable log length range (`minlength`/`maxlength`), and the stem length lost to each cut
(`kerf`):

```@example taper_fit
using DataFrames

products = DataFrame(
    name=["Sawlog", "Pulpwood"],
    sed=[18.0, 8.0],
    minlength=[2.5, 2.0],
    maxlength=[4.0, 3.0],
    kerf=[0.03, 0.03],
)

logassortment(ft, 30.0u"cm", 22.0u"m", products)
```

The greedy algorithm cuts the longest permitted log for the current product until the stem
narrows below its `sed`, then moves to the next product for the remainder of the stem —
`volume`/`logs` report each product's totals as tuples, in table order, alongside the
overall `totalvolume`/`totallogs`.
