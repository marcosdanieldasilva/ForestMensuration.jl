@testset "logassortment" begin
  dbh = repeat([20.0, 25.0, 30.0, 35.0, 40.0, 45.0], inner=4)
  ht = repeat([18.0, 20.0, 22.0, 24.0, 26.0, 28.0], inner=4)
  hr = repeat([0.02, 0.3, 0.6, 0.9], outer=6)
  hi = hr .* ht
  drtrue(hr) = 1.05 - 0.9hr - 0.15hr^2
  di = dbh .* (drtrue.(hr) .+ [0.01sin(3.7i) for i in 1:24])

  ft = fit(Kozak1969(), dbh, ht, hi, di)

  products = DataFrame(
    name=["Sawlog", "Pulpwood"], sed=[18.0, 8.0], minlength=[2.5, 2.0], maxlength=[4.0, 3.0], kerf=[0.03, 0.03],
  )

  result = logassortment(ft, 35.0u"cm", 24.0u"m", products)
  @test result isa DataFrame
  @test nrow(result) == 1
  @test result.product[1] == ("Sawlog", "Pulpwood")
  @test collect(result.volume[1]) ≈ [0.782897, 0.107038] atol = 1e-5
  @test result.logs[1] == (3, 2)
  @test unit(result.totalvolume[1]) == u"m^3"
  @test ustrip(result.totalvolume[1]) ≈ 0.8899350245771176 atol = 1e-6
  @test result.totallogs[1] == 5

  @testset "internal consistency: sums match" begin
    @test sum(result.volume[1]) ≈ ustrip(u"m^3", result.totalvolume[1]) atol = 1e-10
    @test sum(result.logs[1]) == result.totallogs[1]
  end

  @testset "Real and Unitful agree" begin
    resultReal = logassortment(ft, 35.0, 24.0, products)
    @test resultReal.volume[1] == result.volume[1]
    @test resultReal.totalvolume[1] == result.totalvolume[1]
  end

  @testset "removeunits leaves the per-product Tuple columns untouched (already unitless) and strips totalvolume" begin
    stripped = removeunits(result)
    @test "totalvolume (m^3)" ∈ names(stripped)
    @test stripped[1, "totalvolume (m^3)"] ≈ 0.8899350245771176 atol = 1e-6
    @test stripped.volume[1] == result.volume[1]   # already plain numbers, nothing to strip
    @test eltype(stripped.volume) <: Tuple{Vararg{Float64}}
  end

  @testset "kerf reduces usable volume (guards against timbeR's kerf/max-length column mixup)" begin
    zeroKerf = DataFrame(name=["Sawlog", "Pulpwood"], sed=[18.0, 8.0], minlength=[2.5, 2.0], maxlength=[4.0, 3.0], kerf=[0.0, 0.0])
    resultNoKerf = logassortment(ft, 35.0u"cm", 24.0u"m", zeroKerf)
    @test resultNoKerf.totalvolume[1] >= result.totalvolume[1]
    @test resultNoKerf.logs[1] == result.logs[1]   # same length/SED limits, kerf only shifts positions within them
  end

  @testset "stumpheight leaves less stem to buck" begin
    resultStump = logassortment(ft, 35.0u"cm", 24.0u"m", products; stumpheight=0.5u"m")
    @test resultStump.totalvolume[1] < result.totalvolume[1]

    resultStumpReal = logassortment(ft, 35.0, 24.0, products; stumpheight=0.5)
    @test resultStumpReal.totalvolume[1] == resultStump.totalvolume[1]
  end

  @testset "generic to N products (4 here, not just 2)" begin
    products4 = DataFrame(
      name=["Veneer", "Sawlog", "Pulpwood", "Energy"], sed=[25.0, 18.0, 8.0, 4.0],
      minlength=[2.5, 2.5, 2.0, 1.0], maxlength=[3.0, 4.0, 3.0, 2.0], kerf=[0.03, 0.03, 0.03, 0.02],
    )
    result4 = logassortment(ft, 40.0, 26.0, products4)
    @test length(result4.product[1]) == 4
    @test result4.logs[1] == (3, 2, 2, 1)
    @test sum(result4.volume[1]) ≈ ustrip(u"m^3", result4.totalvolume[1]) atol = 1e-10
  end

  @testset "validation" begin
    badProducts = DataFrame(name=["Sawlog"], sed=[18.0])   # missing minlength/maxlength/kerf
    @test_throws ArgumentError logassortment(ft, 35.0, 24.0, badProducts)
    @test_throws DomainError logassortment(ft, 35.0, 24.0, products; stumpheight=30.0)   # >= height
  end
end
