@testset "simplecasualsampling" begin
  v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]

  report = sampling(SimpleCasualSampling, v, 0.05, 10; e=10, α=0.95)
  @test report isa DataFrame
  @test nrow(report) == 1

  @test ustrip(report.vm[1]) ≈ 441.691 atol = 1e-3
  @test report.cv[1] ≈ 10.03 atol = 1e-2
  @test ustrip(report.s2m[1]) ≈ 168.498 atol = 1e-3
  @test ustrip(report.se[1]) ≈ 12.9807 atol = 1e-4
  @test ustrip(report.abserr[1]) ≈ 28.9227 atol = 1e-4
  @test report.relerr[1] ≈ 6.55 atol = 1e-2
  @test ustrip(report.vha[1]) ≈ 8833.82 atol = 1e-2
  @test ustrip(report.vtotal[1]) ≈ 88338.2 atol = 1e-1
  @test ustrip(report.cilower[1]) ≈ 82553.6 atol = 1e-1
  @test ustrip(report.ciupper[1]) ≈ 94122.7 atol = 1e-1
  @test report.pop[1] == "finite"
  @test report.f[1] ≈ 0.945 atol = 1e-3
  @test report.n[1] == 11
  @test report.nreq[1] == 7
  @test report.nmiss[1] == 0
  @test report.N[1] == 200

  # every column that carries a physical unit is a real Unitful quantity column, not a
  # plain number paired with a separate units column -- this is what makes it compatible
  # with `removeunits`/`restoreunits` (ForestFoundations.jl)
  @test eltype(report.vm) <: Unitful.Quantity
  @test eltype(report.vtotal) <: Unitful.Quantity
  @test eltype(report.cv) <: Real
  stripped = removeunits(report)
  @test "vm (m^3)" ∈ names(stripped)
  @test stripped[1, "vm (m^3)"] ≈ 441.691 atol = 1e-3
  restored = restoreunits(stripped)
  @test restored.vm == report.vm

  @testset "units" begin
    reportU = sampling(SimpleCasualSampling, v * u"m^3", 0.05u"ha", 10u"ha")
    @test reportU == report

    reportFt = sampling(SimpleCasualSampling, v * u"ft^3", 0.05u"ac", 10u"ac")
    @test unit(reportFt.vm[1]) == u"ft^3"
    @test unit(reportFt.vha[1]) == u"ft^3/ac"
  end

  @testset "infinite population via 1-reference-area plot" begin
    reportInf = sampling(SimpleCasualSampling, v, 1.0, 1000)
    @test reportInf.pop[1] == "infinite"
  end

  @test_throws ArgumentError sampling(SimpleCasualSampling, [1.0], 0.05, 10)
end
