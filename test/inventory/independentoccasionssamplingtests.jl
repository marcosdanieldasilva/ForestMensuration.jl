@testset "independentoccasionssampling" begin
  v1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0]
  v2 = [24.1, 23.8, 22.9, 25.6, 24.0]

  report = independentoccasionssampling(v1, v2, 0.05, 200, 200)
  @test report isa SamplingReport
  change = report.change

  @test ustrip(change.gm[1]) ≈ 3.91333 atol = 1e-4
  @test ustrip(change.se[1]) ≈ 0.763834 atol = 1e-4
  @test ustrip(change.gtotal[1]) ≈ 782.667 atol = 1e-2

  @test ustrip(report.occasion1.vm[1]) ≈ mean(v1) atol = 1e-9
  @test ustrip(report.occasion2.vm[1]) ≈ mean(v2) atol = 1e-9
  @test report.occasion1.n[1] == 6
  @test report.occasion2.n[1] == 5

  @testset "units" begin
    reportU = independentoccasionssampling(v1 * u"m^3", v2 * u"m^3", 0.05u"ha", 200, 200)
    @test reportU.change == report.change
  end
end
