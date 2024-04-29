module ForestMensuration
  using DataFrames, LinearAlgebra

  include("cubage.jl")

  export cubage, Smalian, Huber, Newton

end
