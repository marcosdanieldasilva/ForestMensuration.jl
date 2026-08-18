@testset "siteclassification.jl" begin
  AGE = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0]
  HDOM = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2, 22.0, 23.1, 24.0, 25.0]

  ageDataPlain() = DataFrame(idade=AGE, hdom=HDOM)
  ageDataUnitful() = DataFrame(idade=AGE, hdom=HDOM .* u"m")

  @testset "siteClassification / hdomClassification round-trip" begin
    m = fit(AllometricModel, @formula(log(hdom) ~ 1 + idade^-1), ageDataPlain())

    site = siteClassification(m, 10.0)
    @test length(site) == length(AGE)
    @test all(>(0), site)
    @test site ≈ [22.8, 19.3, 17.4, 17.8, 17.3, 18.1, 18.0, 18.1, 19.1, 20.1, 20.4, 21.0, 21.5, 22.1] atol = 1e-6

    back = hdomClassification(m, ageDataPlain(), 10.0, site)
    @test all(isapprox.(back, HDOM, atol=1.0))

    @test_throws DomainError siteClassification(m, -1.0)
    @test_throws DomainError siteClassification(m, 0.0)
    @test_throws DomainError hdomClassification(m, ageDataPlain(), 10.0, -site)
  end

  @testset "siteClassification / hdomClassification: Unitful round-trip" begin
    mUnit = fit(AllometricModel, @formula(log(hdom) ~ 1 + idade^-1), ageDataUnitful())

    site = siteClassification(mUnit, 10.0)
    @test nonmissingtype(eltype(site)) <: Quantity
    @test unit(site[1]) == u"m"

    back = hdomClassification(mUnit, ageDataUnitful(), 10.0, site)
    @test nonmissingtype(eltype(back)) <: Quantity
    @test all(isapprox.(ustrip.(back), HDOM, atol=1.0))

    # hdomClassification also accepts a plain-number site vector for a unitful model
    siteBare = ustrip.(site)
    backFromBare = hdomClassification(mUnit, ageDataUnitful(), 10.0, siteBare)
    @test all(isapprox.(ustrip.(backFromBare), ustrip.(back), atol=1e-6))
  end

  @testset "siteTable: automatic breadth reuses frequencytable's class-breadth convention" begin
    m = fit(AllometricModel, @formula(log(hdom) ~ 1 + idade^-1), ageDataPlain())
    table = siteTable(m, 10.0)
    @test table isa DataFrame
    @test nrow(table) == length(unique(AGE))
    @test all(x -> x isa Real, Matrix(table))  # always plain numeric, even from a Unitful model

    # automatic breadth here is 1.0 (frequencytable's "nice number" rounding of the raw
    # Sturges width), and class centers land at a half-breadth offset (_classcenter's
    # convention) — e.g. "S_17.5", not "S_17.6" (the old, unrounded/on-multiple convention)
    @test names(table)[2] == "S_17.5"
    @test ncol(table) == 8

    mUnit = fit(AllometricModel, @formula(log(hdom) ~ 1 + idade^-1), ageDataUnitful())
    tableUnit = siteTable(mUnit, 10.0)
    @test all(x -> x isa Real, Matrix(tableUnit))

    @test_throws DomainError siteTable(m, 10.0, -1.0)
  end

  @testset "siteTable: fixed breadth, same positional convention as frequencytable(x, hi)" begin
    m = fit(AllometricModel, @formula(log(hdom) ~ 1 + idade^-1), ageDataPlain())
    table = siteTable(m, 10.0, 2.0)
    @test table isa DataFrame
    @test names(table) == ["idade", "S_19.0", "S_21.0", "S_23.0"]
  end
end
