## [2.3.0] - 2026-08-19

### Changed

- **Forest inventory sampling extracted to `ForestInventory.jl`:** the 11 sampling designs
  (`sampling`, `SamplingDesign` and its subtypes, `SamplingReport` and its accessors) now
  live in [ForestInventory.jl](https://github.com/JuliaForests/ForestInventory.jl), a new
  package in the `JuliaForests` org. `ForestMensuration.jl` depends on it and re-exports it
  exactly like it already does for `ForestModeling.jl`, so the public API is unchanged —
  `using ForestMensuration` still gives you `sampling(SimpleCasualSampling, ...)` and every
  other inventory function unchanged. `Distributions` and `LinearAlgebra` were dropped from
  `ForestMensuration.jl`'s own dependencies (they were only used by the inventory module).

## [2.0.0] - 2025-02-09

### Breaking Changes

- **Regression Implementation Refactored:**  
  The underlying implementation of the `regression` function has been completely refactored.  
  Previously, it relied on the `GLM.lm` function and its associated structure.  
  It now uses a custom-built regression structure and implementation.  
  **Note:** Although the public API (function call and usage) remains unchanged, this internal change may affect any code that depended on the internal structure of the regression result.

### Other Changes

- Bumped version from 1.0.0 to 2.0.0 to reflect the breaking change in the regression implementation.
