@testset "taperedvolume" begin
  dbh = repeat([20.0, 25.0, 30.0, 35.0, 40.0, 45.0], inner=4)
  ht = repeat([18.0, 20.0, 22.0, 24.0, 26.0, 28.0], inner=4)
  hr = repeat([0.02, 0.3, 0.6, 0.9], outer=6)
  hi = hr .* ht
  drtrue(hr) = 1.05 - 0.9hr - 0.15hr^2
  di = dbh .* (drtrue.(hr) .+ [0.01sin(3.7i) for i in 1:24])

  ft = fit(Kozak1969(), dbh, ht, hi, di)

  wholevol = taperedvolume(ft, 30.0, 22.0)
  @test wholevol isa Unitful.Quantity
  @test unit(wholevol) == u"m^3"
  @test ustrip(wholevol) ≈ 0.6134210338387808 atol = 1e-8

  @testset "3-arg convenience == explicit (0, height)" begin
    explicitWhole = taperedvolume(ft, 30.0u"cm", 22.0u"m", 0.0u"m", 22.0u"m")
    @test explicitWhole == wholevol
  end

  @testset "partial section" begin
    partVol = taperedvolume(ft, 30.0, 22.0, 2.0, 15.0)
    @test ustrip(partVol) ≈ 0.44675757636678787 atol = 1e-8
    @test partVol < wholevol
  end

  @testset "Real and Unitful agree" begin
    vU = taperedvolume(ft, 30.0u"cm", 22.0u"m", 2.0u"m", 15.0u"m")
    @test vU ≈ taperedvolume(ft, 30.0, 22.0, 2.0, 15.0)
  end

  @testset "cross-check against cubage() on a densely-sampled profile" begin
    # different numerical methods (adaptive quadrature on the fitted curve vs.
    # Smalian section sums on discretized samples of the same curve) — agree
    # to within a fraction of a percent, not exactly.
    hs = collect(0.02:0.2:21.98)
    ds = ustrip.(taperdiameter(ft, 30.0, 22.0, hs))
    cb = cubage(hs, ds)
    @test abs(ustrip(cb.vt[1]) - ustrip(wholevol)) / ustrip(wholevol) < 1e-3
  end

  @testset "invalid ranges throw DomainError" begin
    @test_throws DomainError taperedvolume(ft, 30.0, 22.0, 15.0, 2.0)    # hmin > hmax
    @test_throws DomainError taperedvolume(ft, 30.0, 22.0, -1.0, 10.0)   # hmin < 0
    @test_throws DomainError taperedvolume(ft, 30.0, 22.0, 2.0, 30.0)    # hmax > height
  end

  @testset "removeunits / restoreunits" begin
    report = DataFrame(v=wholevol)
    stripped = removeunits(report)
    @test "v (m^3)" ∈ names(stripped)
    @test stripped[1, "v (m^3)"] ≈ ustrip(wholevol) atol = 1e-8
    @test restoreunits(stripped).v == report.v
  end
end
