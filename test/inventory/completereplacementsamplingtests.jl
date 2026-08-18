@testset "completereplacementsampling" begin
  v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 24.6, 17.3]
  v2 = [22.1, 25.8, 23.4, 21.0, 26.6, 20.9, 26.0, 22.5]

  report = sampling(CompleteReplacementSampling, v1, v2, 0.05, 200, 200)
  @test report isa SamplingReport
  chg = change(report)

  @test ustrip(chg.gm[1]) ≈ mean(v2) - mean(v1) atol = 1e-9
  @test ustrip(chg.se[1]) ≈ 0.500612 atol = 1e-4

  @test_throws DimensionMismatch sampling(CompleteReplacementSampling, v1, [1.0, 2.0], 0.05, 200, 200)

  @testset "units" begin
    reportU = sampling(CompleteReplacementSampling, v1 * u"m^3", v2 * u"m^3", 0.05u"ha", 200, 200)
    @test change(reportU) == chg
  end
end
