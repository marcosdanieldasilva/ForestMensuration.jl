@testset "doublesampling" begin
  volume1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, 17.6, 20.8, 21.9]
  volume2 = [22.1, 25.8, 23.4, missing, 26.6, missing, 27.5, missing, missing, 26.0]

  report = doublesampling(volume1, volume2, 0.05, 200)
  @test report isa SamplingReport
  o2 = report.occasion2
  change = report.change

  @test ustrip(o2.vm[1]) ≈ 24.4419 atol = 1e-3
  @test o2.b[1] ≈ 1.1148 atol = 1e-3
  @test o2.m[1] == 6
  @test o2.ntemp[1] == 4
  @test ustrip(change.gm[1]) ≈ 4.00186 atol = 1e-3

  @test_throws ArgumentError doublesampling(volume1, [22.1, 25.8, missing, missing, missing, missing, missing, missing, missing, missing], 0.05, 200)
  @test_throws DimensionMismatch doublesampling(volume1, volume2[1:end-1], 0.05, 200)

  @testset "units" begin
    volume2u = [ismissing(x) ? missing : x * u"m^3" for x in volume2]
    reportU = doublesampling(volume1 * u"m^3", volume2u, 0.05u"ha", 200)
    @test reportU.occasion2 == report.occasion2
  end
end
