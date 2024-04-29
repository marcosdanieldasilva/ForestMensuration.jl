# Tutorial

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