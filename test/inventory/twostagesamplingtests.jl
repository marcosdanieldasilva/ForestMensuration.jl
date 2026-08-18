@testset "twostagesampling" begin
  data = DataFrame(
    primary=repeat(1:5, inner=3),
    volume=[18.2, 19.1, 17.8, 22.4, 23.1, 21.9, 15.1, 16.0, 14.8, 27.3, 28.1, 26.9, 19.8, 20.5, 19.1],
  )

  report = twostagesampling(:primary, :volume, 0.02, 40, 6, data)
  @test report isa SamplingReport
  rt = report.result_table

  @test ustrip(rt.vm[1]) ≈ 20.6733 atol = 1e-3
  @test rt.m[1] == 3
  @test rt.M[1] == 6
  @test rt.n[1] == 5
  @test rt.N[1] == 40

  @testset "units" begin
    reportU = twostagesampling(:primary, :volume, 0.02u"ha", 40, 6, data)
    @test reportU.result_table == report.result_table
  end

  @test_throws ArgumentError twostagesampling(:primary, :volume, 0.02, 40, 2, data)
end
