using Test
using Unitful
using DataFrames
using ForestMensuration

@testset "cubage.jl" begin
  @testset "Diameter interpolation" begin
    h_values = [1.0, 1.3, 2.0, 3.0, 4.0, 5.0]u"m"
    d_values = [30.0, 22.5, 20.2, 15.4, 13.2, 10.9]u"cm"

    @test ustrip(uconvert(u"cm", diameterinterpolation(2.5u"m", h_values, d_values))) ≈ 17.8 atol = 1e-4
    @test ustrip(uconvert(u"cm", diameterinterpolation(3.0u"m", h_values, d_values))) ≈ 15.4 atol = 1e-4
    @test_throws DomainError diameterinterpolation(0.5u"m", h_values, d_values)
  end

  @testset "Height interpolation" begin
    h_values = [1.0, 1.3, 2.0, 3.0, 4.0, 5.0]u"m"
    d_values = [30.0, 22.5, 20.2, 15.4, 13.2, 10.9]u"cm"

    interpolated_height, index = heightinterpolation(17.8u"cm", h_values, d_values)
    @test ustrip(uconvert(u"m", interpolated_height)) ≈ 2.5 atol = 1e-4
    @test index == 4
    @test_throws DomainError heightinterpolation(35.0u"cm", h_values, d_values)
  end

  @testset "Volume helpers" begin
    @test ustrip(uconvert(u"m^3", cylindervolume(18.5u"m", 30.0u"cm"))) ≈ 1.3076879420567515 atol = 1e-6
    @test ustrip(uconvert(u"m^3", conevolume(18.5u"m", 30.0u"cm"))) ≈ 0.4358959806855838 atol = 1e-6
    @test barkfactor([30.0, 22.5, 20.2, 15.4, 13.2, 10.9]u"cm", [1.2, 1.1, 0.85, 0.66, 0.48, 0.0]u"cm") ≈ 0.961764705882353 atol = 1e-6
  end

  @testset "Section methods" begin
    @test ustrip(uconvert(u"m^3", smalian(3.0u"m", 30.0u"cm", 25.0u"cm"))) ≈ 0.17965982987716633 atol = 1e-6
    @test ustrip(uconvert(u"m^3", huber(3.0u"m", 27.5u"cm"))) ≈ 0.1781872083207961 atol = 1e-6
    @test ustrip(uconvert(u"m^3", newton(3.0u"m", 30.0u"cm", 27.5u"cm", 25.0u"cm"))) ≈ 0.1786780821729195 atol = 1e-6
    @test_throws DomainError smalian(-3.0u"m", 30.0u"cm", 25.0u"cm")
  end

  @testset "Form factors and cubage" begin
    vt = 0.3378u"m^3"
    ht = 18.5u"m"
    dbh = 22.7u"cm"
    d_values = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0]u"cm"
    h_values = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8]u"m"

    @test artificialformfactor(vt, ht, dbh) ≈ 0.451176344374475 atol = 1e-6
    @test naturalformfactor(vt, ht, h_values, d_values) ≈ 5.225722786868707 atol = 1e-6
    @test quotientform(ht, dbh, h_values, d_values) ≈ 0.08579295154185025 atol = 1e-6

    result_df = cubage(h_values, d_values; dlimit=2.8u"cm")
    @test isa(result_df, DataFrame)
    @test ustrip(uconvert(u"m^3", result_df.vt[1])) > 0
  end
end