@testset "independentoccasionssampling" begin
  v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0]
  v2 = [24.1, 23.8, 22.9, 25.6, 24.0]

  report = sampling(IndependentOccasionsSampling, v1, v2, 0.05, 200, 200)
  @test report isa SamplingReport
  chg = change(report)

  @test ustrip(chg.gm[1]) ≈ 3.91333 atol = 1e-4
  @test ustrip(chg.se[1]) ≈ 0.763834 atol = 1e-4
  @test ustrip(chg.gtotal[1]) ≈ 782.667 atol = 1e-2

  @test ustrip(occasion1(report).vm[1]) ≈ mean(v1) atol = 1e-9
  @test ustrip(occasion2(report).vm[1]) ≈ mean(v2) atol = 1e-9
  @test occasion1(report).n[1] == 6
  @test occasion2(report).n[1] == 5

  @testset "units" begin
    reportU = sampling(IndependentOccasionsSampling, v1 * u"m^3", v2 * u"m^3", 0.05u"ha", 200, 200)
    @test change(reportU) == chg
  end
end
