## [2.0.0] - 2025-02-09

### Breaking Changes

- **Regression Implementation Refactored:**  
  The underlying implementation of the `regression` function has been completely refactored.  
  Previously, it relied on the `GLM.lm` function and its associated structure.  
  It now uses a custom-built regression structure and implementation.  
  **Note:** Although the public API (function call and usage) remains unchanged, this internal change may affect any code that depended on the internal structure of the regression result.

### Other Changes

- Bumped version from 1.0.0 to 2.0.0 to reflect the breaking change in the regression implementation.
