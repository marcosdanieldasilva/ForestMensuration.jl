using ForestMensuration
using Test
using DataFrames
# Define shared values for tests
@testset "ForestMensuration.jl" begin
  include("cubage_tests.jl")
  include("statistics_tests.jl")
  include("regression_tests.jl")
end
