# ForestMensuration Package:

Documentation for [ForestMensuration](https://github.com/JuliaForests/ForestMensuration.jl).

**ForestMensuration.jl** is a Julia package that offers a comprehensive suite of functions for dendrometric calculations. With a focus on simplicity and efficiency, it streamlines complex forestry computations through an intuitive interface. Key features include:

- **Regression Analysis**: Seamlessly fit linear models—complete with variable transformations and handling of categorical data—to uncover the best relationships in your forestry datasets, including range-safe grouped/stratified variants. Provided by [ForestModeling.jl](https://github.com/JuliaForests/ForestModeling.jl), re-exported here.
- **Tree and Stand Volume Estimation**: Accurately compute tree and stand volumes using a range of methods, including Huber, Smalian, and Newton techniques.
- **Forest Inventory Sampling**: All 11 classic sampling designs — simple, stratified, systematic, cluster/two-stage, horizontal point/Bitterlich, and successive-occasions — each generic to any number of strata/clusters/plots.
- **Stem Taper Equations**: Fit 10 classic published taper forms (Kozak, Bi, Demaerschalk, Max & Burkhart, and more — provided by ForestModeling.jl, re-exported here), then evaluate diameter/height along the stem, integrate volume, or simulate log assortment (sortimentos).
- **Site Productivity Classification**: Guide-curve (anamorphic, delta-method) site index classification from any allometric model fitted with age as its sole continuous regressor.
- **Dendrometric Averaging**: Easily calculate essential metrics such as mean diameter, quadratic mean diameter, and other averages to analyze stand structure.
- **Frequency and Diametric Tables**: Generate detailed frequency and diameter distribution tables to support comprehensive data analysis.

Whether you are conducting research or managing forestry operations, **ForestMensuration.jl** simplifies the process of analyzing dendrometric data, enabling you to perform sophisticated calculations with minimal effort.

```@meta
CurrentModule = ForestMensuration
```

## Installation

The ForestMensuration package is available through the Julia package system and can be added by running Pkg.add("ForestMensuration") or by directly downloading it from the GitHub page:

```julia-repl-repl
pkg> add https://github.com/JuliaForests/ForestMensuration.jl
```

## About the Author

I am Marcos Daniel da Silva, a Forest Engineer with a strong interest in forest mensuration and data analysis. I developed the ForestMensuration.jl package as part of my final course project in Forest Engineering at the Federal University of Santa Maria campus Frederico Westphalen (UFSM), under the guidance of Prof. Dr. Rafaelo Balbinot (UFSM) and Prof. Dr. Alexandre Behling (UFPR). This package represents my dedication to providing accessible tools for forestry professionals and researchers, aiming to simplify complex dendrometric calculations and forest inventory analyses using the Julia programming language.

For more details, please feel free to connect with me on [LinkedIn](https://www.linkedin.com/in/marcosdanieldasilva/?locale=en_US) or by e-mail: [marcosdasilva@5a.tec.br](mailto:marcosdasilva@5a.tec.br).

This package is supported and encouraged by [5A Inteligência e Engenharia](https://5a.tec.br/).

You can also access my final course work through the UFSM repository: [Forestmensuration.jl: Uma Introdução à Aplicações em Julia](https://repositorio.ufsm.br/handle/1/31917?show=full).
