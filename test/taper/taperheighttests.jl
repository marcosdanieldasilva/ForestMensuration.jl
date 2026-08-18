@testset "taperheight" begin
  dbh = repeat([20.0, 25.0, 30.0, 35.0, 40.0, 45.0], inner=4)
  ht = repeat([18.0, 20.0, 22.0, 24.0, 26.0, 28.0], inner=4)
  hr = repeat([0.02, 0.3, 0.6, 0.9], outer=6)
  hi = hr .* ht
  drtrue(hr) = 1.05 - 0.9hr - 0.15hr^2
  di = dbh .* (drtrue.(hr) .+ [0.01sin(3.7i) for i in 1:24])

  ft = fit(Kozak1969(), dbh, ht, hi, di)

  h = taperheight(ft, 30.0, 22.0, 15.0)
  @test h isa Unitful.Quantity
  @test unit(h) == u"m"
  @test ustrip(h) ≈ 12.294953720241237 atol = 1e-6

  @testset "inverts taperdiameter" begin
    dBack = taperdiameter(ft, 30.0, 22.0, ustrip(h))
    @test ustrip(dBack) ≈ 15.0 atol = 1e-6
  end

  @testset "Real and Unitful agree" begin
    hU = taperheight(ft, 30.0u"cm", 22.0u"m", 15.0u"cm")
    @test hU == h
  end

  @testset "unattainable diameter throws DomainError" begin
    # base ≈ 31.5 cm, tip ≈ 0.03 cm for this tree — outside [tip, base] has no solution
    @test_throws DomainError taperheight(ft, 30.0, 22.0, 100.0)   # wider than the base
    @test_throws DomainError taperheight(ft, 30.0, 22.0, 0.001)   # narrower than the tip
  end
end
