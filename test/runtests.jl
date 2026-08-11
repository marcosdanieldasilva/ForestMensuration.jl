using ForestMensuration
using Test
using DataFrames

using ForestMensuration: dm, dg, dw, dz, dd, dh, hm, hd, hg, dominantTreeCount
# Define shared values for tests
@testset "ForestMensuration.jl" begin
  include("cubagetests.jl")
  include("statisticstests.jl")
end
