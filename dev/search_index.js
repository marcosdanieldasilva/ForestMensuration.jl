var documenterSearchIndex = {"docs":
[{"location":"tutorials/#Tutorial","page":"Tutorials","title":"Tutorial","text":"","category":"section"},{"location":"tutorials/#Computing-the-[cubage](@ref)","page":"Tutorials","title":"Computing the cubage","text":"","category":"section"},{"location":"tutorials/#Cubing-a-Simple-Tree","page":"Tutorials","title":"Cubing a Simple Tree","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Cubage for a single tree, with diameters measured in centimeters and heights in meters. The diameter at breast height (DBH) and total height (Ht) must be provided as part of the input values. The CubingMethods can be: Smalian, Newton and Huber.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"using ForestMensuration\n\nd = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0]\n\nh =  [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8]\n\ncubage(Smalian, h, d)","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"With the bark thickness value, it is possible to calculate the bark factor and total and commercial volumes without bark. Note: the provided thickness should be the 'single thickness' in millimeters, i.e., the actual value collected from the field. The function will convert it into 'double thickness'.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"\nbark = [0.9, 0.5, 0.3, 0.2, 0.2, 0.1, 0.0]\n\ncubage(Newton, h, d, bark)","category":"page"},{"location":"tutorials/#Cubing-Multiple-Trees","page":"Tutorials","title":"Cubing Multiple Trees","text":"","category":"section"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"We utilize the function cubage to calculate cubage for multiple trees. The function expects data in a dataframe format where each row represents a tree, and columns contain attributes such as height and diameter. The cubage function can be applied to each tree group using the specified height and diameter columns.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"using ForestMensuration # hide\nusing DataFrames \n\ndata = DataFrame(\n    tree = [148, 148, 148, 148, 148, 148, 148, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222, 222],\n    h = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8, 0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 11.3, 13.3, 15.3, 17.3, 19.5],\n    d = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0, 16.0, 12.0, 11.6, 10.1, 9.4, 8.2, 7.0, 6.0, 4.0, 2.0, 0.0],\n    bark = [0.9, 0.5, 0.3, 0.2, 0.2, 0.1, 0.0, 1.2, 0.5, 0.3, 0.3, 0.2, 0.2, 0.3, 0.0, 0.0, 0.0, 0.0]\n)\n\ncubage(Huber, :tree, :h, :d, data)","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"Additionally, bark thickness values can be provided to calculate bark factors and volumes without bark.","category":"page"},{"location":"tutorials/","page":"Tutorials","title":"Tutorials","text":"cubage(Huber, :tree, :h, :d, :bark, data)","category":"page"},{"location":"reference/#Reference","page":"Reference","title":"Reference","text":"","category":"section"},{"location":"reference/#Contents","page":"Reference","title":"Contents","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/#Index","page":"Reference","title":"Index","text":"","category":"section"},{"location":"reference/","page":"Reference","title":"Reference","text":"Pages = [\"reference.md\"]","category":"page"},{"location":"reference/","page":"Reference","title":"Reference","text":"Modules = [ForestMensuration]","category":"page"},{"location":"reference/#ForestMensuration.Huber","page":"Reference","title":"ForestMensuration.Huber","text":"Huber Method:   The Huber method measures the diameter or circumference at the midpoint of the section, and the volume is determined by:   v = v0 + Σi=1:n(vi) + vt   vi = gi * li   Where:     v0 = volume of the stump;     vi = volume of intermediate sections;     vt = volume of the cone;     g = basal area;     l = length.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ForestMensuration.Newton","page":"Reference","title":"ForestMensuration.Newton","text":"Newton Method:   The Newton method involves measuring at 3 positions along each section (at the ends and in the middle of the logs). Therefore, it is a more laborious method than the others, but the estimated volume will be more accurate.   v = v0 + Σi=1:n(vi) + vt   vi = (gi + gm + gi+1)/2 * li   Where:     v0 = volume of the stump;     vi = volume of intermediate sections;     vt = volume of the cone;     g = basal area;     gm = basal area at the midpoint of the section;     l = length.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ForestMensuration.Smalian","page":"Reference","title":"ForestMensuration.Smalian","text":"Smalian Method:   The Smalian method measures diameters or circumferences at the ends of each section and calculates the total volume by:   Vt = v0 + Σi=1:n(vi) + vt   v0 = g0 * l0   vi = (gi+gi+1)/2 * li   vt = (1/3) * gn * ln   Where:     v0 = volume of the stump;     vi = volume of intermediate sections;     vt = volume of the cone;     g = basal area;     l = length.\n\n\n\n\n\n","category":"type"},{"location":"reference/#ForestMensuration.cubage-Tuple{Type{<:ForestMensuration.CubingMethod}, Vector{<:Real}, Vector{<:Real}}","page":"Reference","title":"ForestMensuration.cubage","text":"Calculate tree cubage using Smalian, Newton, or Huber methods.\n\nThe methods involve dividing the tree trunk into n sections (logs). In each section, diameters and lengths are measured at positions that vary according to the technique employed. Thus, the volume of the sections and the total volume are determined by summing the volume of the sections. Determination can be carried out on felled trees or standing trees using equipment such as the Bitterlich relascope.\n\nForm Factor (ff): The form factor is a reduction factor used to determine the volume of a standing tree. It is the ratio of the rigorous volume of the tree to the volume of a reference cylinder. The volume is calculated by the product of the form factor, basal area, and total height. v = g + h + f Where:   g = basal area;   h = total height;   f = form factor, either natural or artificial.\nArtificial Form Factor (aff: f₁,₃): For the calculation of the artificial form factor, the volume of the reference cylinder will have a diameter equal to the tree's DBH. f1,3 = Rigorous Vol / Cylinder Vol 1,3 Where:   Rigorous Vol = total volume determined by one of the methods: Smalian, Huber, or Newton;   Cylinder Vol 1,3 = volume of a cylinder with height and diameter equal to the total height and DBH of the tree.\nNatural Form Factor (nff: f₀,₁ₕ): For the calculation of the natural form factor, the volume of the reference cylinder will have a diameter equal to the diameter taken at 1/10 of the total height. f0,1h = Rigorous Vol / Cylinder Vol 0,1 Where:   Rigorous Vol = total volume determined by one of the methods: Smalian, Huber, or Newton;   Cylinder Vol 0,1 = volume of a cylinder with height equal to the total height of the tree and diameter taken at 1/10 of the total height.   Interpolate diameter at a given height using linear interpolation.\nForm Quotient (qf: Q):   The natural decrease in diameter along the trunk defines the so-called form quotient, which is a ratio between diameters. An example of a form quotient is the Schiffel form quotient, given by:\nQ = D(1/2H) / DBH, Where     Q < 1     D(1/2H) = diameter measured at half the total height of the tree.\nSimilar to the form factor, the volume of a tree, with or without bark, can be obtained by multiplying the volume of a cylinder by the average form quotient, suitable for the species and the desired volume to be estimated.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = ForestMensuration","category":"page"},{"location":"#ForestMensuration","page":"Home","title":"ForestMensuration","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ForestMensuration.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"#Install","page":"Home","title":"Install","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"\"pkg> add https://github.com/marcosdanieldasilva/ForestMensuration.jl\"","category":"page"}]
}
