@testset "systematicsampling" begin
  v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]

  report = systematicsampling(v, 0.05, 10)
  @test nrow(report) == 1

  @test ustrip(report.vm[1]) ≈ 441.691 atol = 1e-3
  @test ustrip(report.s2m[1]) ≈ 108.254 atol = 1e-2
  @test ustrip(report.se[1]) ≈ 10.4045 atol = 1e-3
  @test report.k[1] == 1
  @test report.n[1] == 11
  @test report.N[1] == 200

  @testset "with lines" begin
    line = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2]
    reportLine = systematicsampling(v, 0.05, 10; line=line)
    @test reportLine.k[1] == 2
    # the difference straddling the two lines is excluded, so the variance differs from
    # the single-line case
    @test ustrip(reportLine.s2m[1]) != ustrip(report.s2m[1])
  end

  @testset "units" begin
    reportU = systematicsampling(v * u"m^3", 0.05u"ha", 10u"ha")
    @test reportU == report
  end

  @test_throws DimensionMismatch systematicsampling(v, 0.05, 10; line=[1, 2, 3])
  @test_throws ArgumentError systematicsampling([1.0], 0.05, 10)
end
