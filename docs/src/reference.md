```@meta
CurrentModule = ForestMensuration
```

Documentation for [ForestMensuration](https://github.com/JuliaForests/ForestMensuration.jl).

# Reference

```@docs
ForestMensuration
ForestInventory
ForestModeling
```

## Cubage

```@docs
cubage
smalian
newton
huber
hohenadl
artificialformfactor
naturalformfactor
quotientform
barkfactor
conevolume
cylindervolume
diameterinterpolation
heightinterpolation
```

## Dendrometric Averages

```@docs
dmetrics
hmetrics
standmetrics
dm
dg
dw
dz
dd
dh
hm
hd
hg
```

## Frequency and Diametric Tables

```@docs
frequencytable
diametrictable
```

## Forest Inventory Sampling (re-exported from ForestInventory.jl)

```@autodocs
Modules = [ForestInventory]
Order = [:function, :type, :constant]
```

## Site Classification

```@docs
siteClassification
hdomClassification
siteTable
```

## Stem Taper Application

```@docs
taperdiameter
taperheight
taperedvolume
logassortment
```

## Regression and Stem Taper Fitting (re-exported from ForestModeling.jl)

```@autodocs
Modules = [ForestModeling]
Order = [:function, :type, :constant]
```
