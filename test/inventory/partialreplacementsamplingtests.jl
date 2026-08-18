@testset "partialreplacementsampling" begin
  volume1 = [18.2, 21.4, 19.8, 20.1, 22.5, 19.0, 23.1, missing, missing]
  volume2 = [missing, missing, 23.4, 24.0, 26.6, 22.9, 27.5, 25.2, 21.8]

  report = partialreplacementsampling(volume1, volume2, 0.05, 200)
  @test report isa SamplingReport
  o2 = occasion2(report)
  chg = change(report)

  @test ustrip(o2.vm[1]) ≈ 24.3459 atol = 1e-3
  @test o2.m[1] == 5
  @test o2.u[1] == 2
  @test o2.v[1] == 2
  @test 0 <= o2.c[1] <= 1
  @test ustrip(chg.gm[1]) ≈ 3.76019 atol = 1e-3

  @test_throws DimensionMismatch partialreplacementsampling(volume1, volume2[1:end-1], 0.05, 200)
  @test_throws ArgumentError partialreplacementsampling(
    [1.0, 2.0, missing, missing], [missing, missing, 3.0, 4.0], 0.05, 200)

  @testset "units" begin
    volume1u = [ismissing(x) ? missing : x * u"m^3" for x in volume1]
    volume2u = [ismissing(x) ? missing : x * u"m^3" for x in volume2]
    reportU = partialreplacementsampling(volume1u, volume2u, 0.05u"ha", 200)
    @test occasion2(reportU) == o2
  end
end
