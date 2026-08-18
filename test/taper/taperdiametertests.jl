@testset "taperdiameter" begin
  dbh = repeat([20.0, 25.0, 30.0, 35.0, 40.0, 45.0], inner=4)
  ht = repeat([18.0, 20.0, 22.0, 24.0, 26.0, 28.0], inner=4)
  hr = repeat([0.02, 0.3, 0.6, 0.9], outer=6)
  hi = hr .* ht
  drtrue(hr) = 1.05 - 0.9hr - 0.15hr^2
  di = dbh .* (drtrue.(hr) .+ [0.01sin(3.7i) for i in 1:24])

  ft = fit(Kozak1969(), dbh, ht, hi, di)

  d = taperdiameter(ft, 30.0, 22.0, 5.0)
  @test d isa Unitful.Quantity
  @test unit(d) == u"cm"
  @test ustrip(d) ≈ 25.127206142774387 atol = 1e-8

  @testset "Real and Unitful agree" begin
    dU = taperdiameter(ft, 30.0u"cm", 22.0u"m", 5.0u"m")
    @test dU == d
    dMixed = taperdiameter(ft, 300.0u"mm", 22.0u"m", 5.0u"m")   # 300 mm == 30 cm
    @test ustrip(dMixed) ≈ ustrip(d) atol = 1e-8
  end

  @testset "vector h" begin
    dvec = taperdiameter(ft, 30.0, 22.0, [1.0, 5.0, 10.0])
    @test length(dvec) == 3
    @test dvec[2] == d
    @test ustrip.(dvec) ≈ [30.269996225036028, 25.127206142774387, 18.290148546680445] atol = 1e-8
    @test issorted(ustrip.(dvec), rev=true)   # diameter decreases with height over this range

    dvecU = taperdiameter(ft, 30.0u"cm", 22.0u"m", [1.0, 5.0, 10.0] .* u"m")
    @test dvecU == dvec
  end
end
