using ForestMensuration
using Test
using DataFrames
using Statistics
using Unitful

using ForestMensuration: dm, dg, dw, dz, dd, dh, hm, hd, hg, dominantTreeCount
# Define shared values for tests
@testset "ForestMensuration.jl" begin
  include("cubagetests.jl")
  include("statisticstests.jl")
  include("tapertests.jl")
  include("siteclassificationtests.jl")
end
