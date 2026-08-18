@testset "multistartsystematicsampling" begin
  data = DataFrame(
    start=repeat(1:4, inner=5),
    volume=[18.2, 19.1, 17.8, 18.9, 19.4, 22.4, 23.1, 21.9, 22.8, 23.4,
      15.1, 16.0, 14.8, 15.6, 15.9, 20.5, 21.2, 19.8, 20.6, 21.0],
  )

  report = multistartsystematicsampling(:start, :volume, 0.02, 15, data)
  @test report isa SamplingReport
  rt = resultTable(report)

  @test ustrip(rt.vm[1]) ≈ 19.375 atol = 1e-3
  @test rt.M[1] == 5
  @test rt.n[1] == 4
  @test rt.N[1] == 150

  @testset "units" begin
    reportU = multistartsystematicsampling(:start, :volume, 0.02u"ha", 15u"ha", data)
    @test resultTable(reportU) == rt
  end
end
