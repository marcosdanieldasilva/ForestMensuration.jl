@testset "stratifiedsampling" begin
  data = DataFrame(
    stratum=[1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3],
    volume=[18.2, 21.4, 19.8, 20.1, 32.5, 35.1, 30.8, 12.4, 11.9, 13.6, 12.8, 13.1],
  )

  report = stratifiedsampling(:stratum, :volume, 0.1, [12.0, 8.0, 20.0], data)
  @test report isa SamplingReport
  rt = resultTable(report)
  @test nrow(rt) == 1

  @test ustrip(rt.vm[1]) ≈ 18.9025 atol = 1e-3
  @test ustrip(rt.vtotal[1]) ≈ 7561.0 atol = 1e-1
  @test rt.nh[1] == (4, 3, 5)
  @test rt.n[1] == 12
  @test rt.N[1] == 400
  @test eltype(rt.vm) <: Unitful.Quantity

  @test nrow(auxiliaryTable(report)) == 3
  @test auxiliaryTable(report).n == [4, 3, 5]
  @test ustrip(auxiliaryTable(report).x̅) ≈ [19.875, 32.8, 12.76] atol = 1e-3

  @test nrow(anova(report)) == 3
  @test anova(report).DOF == [2, 9, 11]

  @testset "one stratum matches simplecasualsampling" begin
    v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]
    data1 = DataFrame(stratum=fill(1, 11), volume=v)
    r1 = resultTable(stratifiedsampling(:stratum, :volume, 0.05, [10.0], data1))
    r2 = simplecasualsampling(v, 0.05, 10)
    @test ustrip(r1.vm[1]) ≈ ustrip(r2.vm[1]) atol = 1e-6
    @test ustrip(r1.s2m[1]) ≈ ustrip(r2.s2m[1]) atol = 1e-6
    @test ustrip(r1.se[1]) ≈ ustrip(r2.se[1]) atol = 1e-6
    @test r1.nreq[1] == r2.nreq[1]
    @test r1.pop[1] == r2.pop[1]
    @test r1.f[1] ≈ r2.f[1] atol = 1e-6
    # the ANOVA between-strata line is undefined with a single stratum (0 numerator DOF)
    r1full = stratifiedsampling(:stratum, :volume, 0.05, [10.0], data1)
    @test ismissing(anova(r1full).F[1])
  end

  @testset "units" begin
    reportU = stratifiedsampling(:stratum, :volume, 0.1u"ha", [12.0, 8.0, 20.0]u"ha", data)
    @test resultTable(reportU) == rt
  end

  @test_throws ArgumentError stratifiedsampling(:stratum, :volume, 0.1, [12.0, 8.0], data)
end
