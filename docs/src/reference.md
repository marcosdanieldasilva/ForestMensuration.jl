```@meta
CurrentModule = ForestMensuration
```

Documentation for [ForestMensuration](https://github.com/JuliaForests/ForestMensuration.jl).

# Reference

```@docs
ForestMensuration
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

## Forest Inventory Sampling

```@docs
SamplingReport
simplecasualsampling
stratifiedsampling
systematicsampling
multistartsystematicsampling
clustersampling
twostagesampling
independentoccasionssampling
completereplacementsampling
partialreplacementsampling
doublesampling
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
