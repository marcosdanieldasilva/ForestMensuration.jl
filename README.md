# ForestMensuration

ForestMensuration.jl provides functions for dendrometric and forest inventory calculations in Julia.

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://marcosdanieldasilva.github.io/ForestMensuration.jl/stable/forestmensuration/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://marcosdanieldasilva.github.io/ForestMensuration.jl/dev/forestmensuration/)
[![Build Status](https://github.com/marcosdanieldasilva/ForestMensuration.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/marcosdanieldasilva/ForestMensuration.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/marcosdanieldasilva/ForestMensuration.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/marcosdanieldasilva/ForestMensuration.jl)

## Installation

```julia
using Pkg
Pkg.add("ForestMensuration")
```

## Example Usage

### Cubage of a Single Tree

Compute the volume of a single tree using the Smalian method.

```julia
using ForestMensuration

# Diameters at different heights (cm)
d = [30.0, 25.0, 20.0, 15.0, 10.0, 5.0, 0.0];

# Corresponding heights (m)
h = [0.7, 1.3, 2.0, 4.0, 6.0, 8.0, 10.0];

# Calculate cubage using the Smalian method
cubage(Smalian, h, d)
1×11 DataFrame
 Row │ vt        v0         vc        vr       vn        dbh      ht       hc       aff       nff       qf
     │ Float64   Float64    Float64   Float64  Float64   Float64  Float64  Float64  Float64   Float64   Float64
─────┼──────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ 0.199328  0.0494801  0.148538      0.0  0.001309     25.0     10.0      8.0  0.406067  0.335592      0.5
```

### Regression Example

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
models = regression(:h, :d, data);

# Select the best model
top_model = criteria_selection(models)
log_minus(h - 1.3) = 4.1496 - 37.62 * d ^ -1
```
