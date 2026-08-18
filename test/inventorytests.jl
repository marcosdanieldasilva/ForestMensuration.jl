@testset "forest inventory sampling" begin
  include("inventory/simplecasualsamplingtests.jl")
  include("inventory/stratifiedsamplingtests.jl")
  include("inventory/systematicsamplingtests.jl")
  include("inventory/clustersamplingtests.jl")
  include("inventory/horizontalpointsamplingtests.jl")
  include("inventory/multistartsystematicsamplingtests.jl")
  include("inventory/twostagesamplingtests.jl")
  include("inventory/independentoccasionssamplingtests.jl")
  include("inventory/completereplacementsamplingtests.jl")
  include("inventory/partialreplacementsamplingtests.jl")
  include("inventory/doublesamplingtests.jl")
end
