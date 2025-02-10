# ForestMensuration

ForestMensuration.jl provides advanced functions for dendrometric calculations in Julia.  
Its focus is on accurate **tree cubage (volume estimation)**, modeling **hipometric relationships** (e.g. diameter versus height and dominant height versus age for site index estimation), as well as **volumetry, allometry,** and **biomass** estimation. These methods are essential for forest mensuration and biometrics, supporting forest inventory, management, and research.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://marcosdanieldasilva.github.io/ForestMensuration.jl/stable/forestmensuration/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://marcosdanieldasilva.github.io/ForestMensuration.jl/dev/forestmensuration/)
[![Build Status](https://github.com/marcosdanieldasilva/ForestMensuration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/marcosdanieldasilva/ForestMensuration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/marcosdanieldasilva/ForestMensuration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/marcosdanieldasilva/ForestMensuration.jl)

## Installation

Install the package via Julia's package manager:

```julia
using Pkg
Pkg.add("ForestMensuration")
```

## Overview

ForestMensuration.jl is designed for professionals and researchers in forestry, dendrometry, and forest biometrics. Its key features include:

- **Cubage (Tree Volume Estimation):**  
  Provides rigorous methods (e.g., the Smalian method) for calculating tree volume.

- **Hipometric Relationships:**  
  Implements regression-based methods for modeling relationships such as:

  - Diameter as a function of height.
  - Dominant height as a function of age for site index calculations.

- **Allometry and Biomass Estimation:**  
  Offers tools to predict tree volume and biomass based on measurements like diameter at breast height and total height.

These functions are crucial for generating accurate forest inventories and supporting decision-making in forest management.

## Example Usage

### Cubage Calculation

Compute the volume of a single tree using the Smalian method.

```julia
using ForestMensuration

# Diameters at different heights (in cm)
d = [30.0, 25.0, 20.0, 15.0, 10.0, 5.0, 0.0]

# Corresponding heights (in m)
h = [0.7, 1.3, 2.0, 4.0, 6.0, 8.0, 10.0]

# Calculate tree cubage using the Smalian method
cubage(Smalian, h, d)
1×11 DataFrame
 Row │ vt        v0         vc        vr       vn        d        h        hc       aff       nff       qf
     │ Float64   Float64    Float64   Float64  Float64   Float64  Float64  Float64  Float64   Float64   Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 0.199328  0.0494801  0.148538      0.0  0.001309     25.0     10.0      8.0  0.406067  0.335592      0.5
```

### Regression Analysis for Hipometric Relationships

Fit a regression model between tree height (`h`) and diameter (`d`).

```julia
using ForestMensuration
using DataFrames

# Sample data with 10 observations
data = DataFrame(
    h = [10.2, 15.3, 14.8, 9.7, 16.5, 13.1, 11.6, 12.4, 14.2, 15.0],
    d = [20.5, 25.3, 24.1, 18.7, 26.2, 22.5, 19.8, 21.0, 23.4, 24.5]
)

# Perform regression analysis between height and diameter
fitted_hipso = regression(:h, :d, data);

# You can evaluate and rank them based on specific criteria using the criteria_table function.
criteria_table(fitted_hipso)
10×9 DataFrame
 Row │ model                              rank   adjr2     d         syx      aic      bic      normality  significance
     │ LinearMo…                          Int64  Float64   Float64   Float64  Float64  Float64  Float64    Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ log_minus(h - 1.3) = 4.1496 - 37…      7  0.935832  0.985189  4.33559  21.3377  21.9429        1.0           1.0
   2 │ log(h) = 4.0815 - 33.72 * d ^ -1      12  0.935775  0.985139  4.3375   21.3465  21.9517        1.0           1.0
   3 │ 1 / √(h - 1.3) = 0.1685 + 61.34 …     17  0.935726  0.985125  4.33915  21.3541  21.9593        1.0           1.0
   4 │ log1p(h) = 4.0449 - 31.24 * d ^ …     22  0.935663  0.985092  4.34128  21.3639  21.9691        1.0           1.0
   5 │ h = -6.5776 + 0.8787 * d              27  0.935641  0.985057  4.34202  21.3674  21.9725        1.0           1.0
   6 │ h = -17.3383 + 3.157 * log(d) ^ 2     32  0.935498  0.985022  4.34685  21.3896  21.9948        1.0           1.0
   7 │ 1 / √h = 0.1721 + 51.94 * d ^ -2      37  0.935405  0.984988  4.34999  21.404   22.0092        1.0           1.0
   8 │ h = -47.6116 + 19.56 * log(d)         44  0.934714  0.984833  4.37318  21.5104  22.1155        1.0           1.0
   9 │ log1p(h) = -1.7544 + 1.414 * log…     46  0.934611  0.984971  4.37663  21.5261  22.1313        1.0           1.0
  10 │ log(h) = -2.1777 + 1.526 * log(d)     51  0.934183  0.984899  4.39093  21.5914  22.1966        1.0           1.0

# Os simply select the best model based on predefined criteria
top_model = criteria_selection(fitted_hipso)
log_minus(h - 1.3) = 4.1496 - 37.62 * d ^ -1
```

### Volume Regression Example

Use volume data to build a regression model. The sample data includes tree volume (`v`), total height (`h`), and diameter at breast height (`d`).

```julia
using ForestMensuration
using DataFrames

# Sample volume data
volume = DataFrame(
  v = [0.153963, 0.151394, 0.178328, 0.221116, 0.225294, 0.183733, 0.210808, 0.200058, 0.183102, 0.231159, 0.272580, 0.224460, 0.322517, 0.246280, 0.264114, 0.308264, 0.326833, 0.316227, 0.331639, 0.354205, 0.421524, 0.331393, 0.421420, 0.462966, 0.427326, 0.482238, 0.345795, 0.476053, 0.556150, 0.441768, 0.543030, 0.449364, 0.450939, 0.520234, 0.422021, 0.410333, 0.598158, 0.763291, 0.610406, 0.596743, 0.667189, 0.653583, 0.626541, 0.638496, 0.683411, 0.740718, 0.710917, 0.737058, 0.714613, 0.749651, 0.786985, 0.781447, 0.727167, 0.960877, 0.768727, 0.690346, 0.971356, 0.989444, 0.809371, 0.945937, 0.892852, 1.074173, 1.176289, 1.003187, 1.010025, 1.054886, 1.030369, 1.176427, 1.021023, 1.098403, 1.050754, 1.070469, 1.118549, 1.559863, 1.083926, 1.322375, 1.271627, 1.057683, 1.262937, 1.496877, 1.398170, 1.317415, 1.653541, 1.313996, 1.611356, 1.399721, 1.300030, 1.480614, 1.583871, 1.454267, 1.616712],
  h = [18.4, 18.9, 20.2, 20.4, 22.5, 19.2, 21.8, 20.0, 19.2, 22.3, 23.0, 19.0, 23.4, 19.1, 20.2, 21.9, 22.2, 21.8, 23.4, 23.2, 22.4, 21.5, 23.5, 23.5, 24.4, 22.9, 21.2, 22.7, 22.6, 22.4, 23.6, 20.7, 22.5, 23.7, 21.0, 22.0, 24.0, 25.6, 23.0, 24.0, 23.3, 23.2, 23.1, 23.0, 23.4, 26.0, 24.1, 23.8, 24.6, 24.2, 24.5, 25.3, 22.1, 25.2, 25.0, 22.6, 24.2, 26.0, 23.5, 25.3, 24.8, 25.1, 25.0, 23.8, 23.0, 24.4, 25.0, 25.1, 24.9, 24.8, 24.2, 23.7, 25.2, 29.2, 23.7, 25.1, 24.9, 24.2, 25.1, 26.7, 25.5, 25.5, 26.2, 23.9, 25.0, 26.5, 26.9, 25.2, 26.0, 26.5, 25.2],
  d = [16.9, 17.2, 17.5, 17.8, 17.8, 17.8, 17.8, 18.5, 18.5, 18.8, 19.1, 19.7, 19.7, 19.7, 20.1, 20.1, 21.0, 21.0, 21.3, 22.0, 22.0, 22.3, 22.6, 22.6, 23.6, 23.6, 23.6, 23.6, 24.2, 24.5, 24.5, 24.5, 25.1, 25.5, 25.5, 26.1, 26.4, 27.1, 27.1, 27.4, 27.4, 27.7, 28.0, 28.0, 28.3, 28.6, 29.3, 29.6, 29.9, 30.2, 30.2, 31.2, 31.2, 31.5, 31.5, 31.5, 32.1, 32.1, 32.5, 33.1, 33.4, 33.7, 33.7, 34.1, 34.4, 35.0, 35.0, 35.0, 35.0, 35.7, 36.3, 36.3, 36.6, 37.2, 37.2, 37.9, 37.9, 38.2, 38.2, 38.8, 38.8, 39.5, 39.8, 39.8, 40.1, 40.1, 40.1, 41.1, 41.4, 41.7, 41.7],
)

# Build a fitted regression models for volume
fitted_volume = regression(:v, :h, :d, volume);

# You can evaluate and rank them based on specific criteria using the criteria_table function.
criteria_table(fitted_volume)
10×9 DataFrame
 Row │ model                              rank   adjr2     d         syx      aic       bic       normality  significance
     │ LinearMo…                          Int64  Float64   Float64   Float64  Float64   Float64   Float64    Float64
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ log(v) = -10.2033 + 0.005545 * h…      7  0.973721  0.99354   9.13789  -219.724  -209.68         1.0           1.0
   2 │ log(v) = -10.4441 + 0.005997 * h…     12  0.973696  0.993532  9.14211  -219.64   -209.596        1.0           1.0
   3 │ log(v) = -6.5575 + 0.495 * d + 0…     18  0.9736    0.993491  9.1589   -219.306  -209.262        1.0           1.0
   4 │ log1p(v) = -0.153 + 0.0008139 * …     23  0.973571  0.993482  9.16383  -219.208  -209.165        1.0           1.0
   5 │ log(v) = -7.4964 + 0.002674 * d …     29  0.973545  0.993471  9.16841  -219.117  -209.073        1.0           1.0
   6 │ log1p(v) = -0.1395 + 0.0004409 *…     32  0.973533  0.993473  9.17053  -219.075  -209.031        1.0           1.0
   7 │ log(v) = -12.5552 + 0.2311 * h +…     33  0.973515  0.993499  9.17355  -219.015  -208.972        1.0           1.0
   8 │ log(v) = -7.4379 + 0.0007277 * d…     42  0.973501  0.993461  9.17605  -218.965  -208.922        1.0           1.0
   9 │ log(v) = -10.1869 + 0.09427 * h …     47  0.973454  0.993452  9.18411  -218.806  -208.762        1.0           1.0
  10 │ log1p(v) = -0.1662 + 0.002532 * …     52  0.97339   0.993441  9.19527  -218.585  -208.541        1.0           1.0
```

### Biomass Regression Example

Fit a regression model using biomass data. The data includes tree biomass (`p`), total height (`h`), and diameter at breast height (`d`).

```julia
using ForestMensuration
using DataFrames

# Sample biomass data
biomass = DataFrame(
  d = [6.5, 8.0, 8.5, 8.5, 9.5, 9.5, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.2, 10.5, 10.5, 10.5, 11.0, 11.0, 11.0, 11.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.5, 12.5, 13.0, 13.0, 13.0, 13.0, 13.0, 13.5, 13.5, 14.0, 14.0, 14.0, 14.5, 14.5, 14.5, 14.5, 14.5, 15.0, 15.0, 15.5, 16.0],
  h = [11.5, 15.2, 14.8, 17.0, 16.9, 19.8, 17.7, 19.0, 18.2, 15.8, 17.5, 17.8, 16.7, 18.0, 17.9, 16.0, 18.8, 17.1, 19.3, 19.5, 18.8, 20.0, 19.1, 16.8, 17.9, 19.7, 19.0, 19.9, 20.5, 20.1, 16.1, 19.3, 20.2, 19.2, 21.1, 21.2, 21.8, 18.9, 20.6, 20.6, 21.8, 22.2, 19.4, 20.5, 21.1, 21.3],
  p = [5.832, 18.825, 16.065, 15.657, 23.476, 29.376, 25.247, 31.962, 31.954, 27.356, 25.163, 29.865, 23.391, 32.887, 38.719, 31.059, 31.689, 26.379, 36.608, 44.182, 37.616, 40.209, 45.496, 39.182, 40.036, 58.883, 56.035, 64.647, 65.173, 59.035, 33.082, 50.726, 53.361, 56.979, 77.001, 84.832, 77.141, 65.280, 65.389, 65.240, 76.359, 69.406, 71.757, 81.424, 100.448, 89.240]
)

# Build a fitted regression models for biomass
fitted_biomass = regression(:p, :h, :d, biomass);

# You can evaluate and rank them based on specific criteria using the criteria_table function.
criteria_table(fitted_biomass)
10×9 DataFrame
 Row │ model                              rank   adjr2     d         syx      aic      bic      normality  significance
     │ LinearMo…                          Int64  Float64   Float64   Float64  Float64  Float64  Float64    Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ log1p(p) = -0.8833 + 0.6421 * lo…     14  0.933133  0.983143  12.3616  296.749  300.406        1.0           1.0
   2 │ log(p) = -0.8808 + 0.1646 * h & …     15  0.934704  0.983424  12.2155  297.655  303.141        1.0           1.0
   3 │ log1p(p) = 1.0968 - 0.4876 * d -…     23  0.93313   0.983627  12.3619  300.751  308.066        1.0           1.0
   4 │ log(p) = 0.9015 - 0.5357 * d - 1…     29  0.932835  0.983545  12.3892  300.954  308.268        1.0           1.0
   5 │ log(p) = -5.6135 + 1.587 * log(h…     30  0.932394  0.983012  12.4297  299.254  304.74         1.0           1.0
   6 │ log(p) = -1.0861 + 0.6663 * log(…     31  0.93192   0.983075  12.4732  297.576  301.233        1.0           1.0
   7 │ p = 137.2894 - 9.591 * h - 53.93…     32  0.932762  0.983547  12.3959  301.004  308.318        1.0           1.0
   8 │ log1p(p) = -5.1686 + 1.48 * log(…     38  0.932089  0.982721  12.4578  299.462  304.948        1.0           1.0
   9 │ p = 90.9551 - 56.86 * d - 4.853 …     40  0.932346  0.983436  12.4342  301.288  308.602        1.0           1.0
  10 │ log1p(p) = 0.8836 + 0.01769 * h …     43  0.93165   0.982827  12.498   299.758  305.244        1.0           1.0
```

## Keywords

Forest mensuration, dendrometry, forest inventory, tree cubage, hipometric relationship, volumetry, allometry, biomass, forest biometrics.

## License

This project is licensed under the MIT License.
