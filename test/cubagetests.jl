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

    # unitless convenience: heights in m, diameters in cm
    @test diameterinterpolation(2.5, [1.0, 1.3, 2.0, 3.0, 4.0, 5.0], [30.0, 22.5, 20.2, 15.4, 13.2, 10.9]) ==
          diameterinterpolation(2.5u"m", h_values, d_values)
  end

  @testset "Height interpolation" begin
    h_values = [1.0, 1.3, 2.0, 3.0, 4.0, 5.0]u"m"
    d_values = [30.0, 22.5, 20.2, 15.4, 13.2, 10.9]u"cm"

    interpolated_height, index = heightinterpolation(17.8u"cm", h_values, d_values)
    @test ustrip(uconvert(u"m", interpolated_height)) ≈ 2.5 atol = 1e-4
    @test index == 4
    @test_throws DomainError heightinterpolation(35.0u"cm", h_values, d_values)

    # unitless convenience: diameters in cm, heights in m
    @test heightinterpolation(17.8, [1.0, 1.3, 2.0, 3.0, 4.0, 5.0], [30.0, 22.5, 20.2, 15.4, 13.2, 10.9]) ==
          heightinterpolation(17.8u"cm", h_values, d_values)
  end

  @testset "Volume helpers" begin
    @test ustrip(uconvert(u"m^3", cylindervolume(18.5u"m", 30.0u"cm"))) ≈ 1.3076879420567515 atol = 1e-6
    @test ustrip(uconvert(u"m^3", conevolume(18.5u"m", 30.0u"cm"))) ≈ 0.4358959806855838 atol = 1e-6
    @test barkfactor([30.0, 22.5, 20.2, 15.4, 13.2, 10.9]u"cm", [1.2, 1.1, 0.85, 0.66, 0.48, 0.0]u"cm") ≈ 0.961764705882353 atol = 1e-6

    # zero is a legitimate stump height / tip diameter, not a domain error
    @test cylindervolume(0.0u"m", 30.0u"cm") == 0.0u"m^3"
    @test conevolume(1.0u"m", 0.0u"cm") == 0.0u"m^3"
    @test_throws DomainError cylindervolume(-1.0u"m", 30.0u"cm")
    @test_throws DomainError conevolume(1.0u"m", -1.0u"cm")

    # unitless convenience: heights/diameters in m/cm
    @test cylindervolume(18.5, 30.0) == cylindervolume(18.5u"m", 30.0u"cm")
    @test conevolume(18.5, 30.0) == conevolume(18.5u"m", 30.0u"cm")
    @test barkfactor([30.0, 22.5], [1.2, 1.1]) == barkfactor([30.0, 22.5]u"cm", [1.2, 1.1]u"cm")
  end

  @testset "Section methods" begin
    @test ustrip(uconvert(u"m^3", smalian(3.0u"m", 30.0u"cm", 25.0u"cm"))) ≈ 0.17965982987716633 atol = 1e-6
    @test ustrip(uconvert(u"m^3", huber(3.0u"m", 27.5u"cm"))) ≈ 0.1781872083207961 atol = 1e-6
    @test ustrip(uconvert(u"m^3", newton(3.0u"m", 30.0u"cm", 27.5u"cm", 25.0u"cm"))) ≈ 0.1786780821729195 atol = 1e-6
    @test_throws DomainError smalian(-3.0u"m", 30.0u"cm", 25.0u"cm")

    # unitless convenience: log length in m, diameters in cm
    @test smalian(3.0, 30.0, 25.0) == smalian(3.0u"m", 30.0u"cm", 25.0u"cm")
    @test huber(3.0, 27.5) == huber(3.0u"m", 27.5u"cm")
    @test newton(3.0, 30.0, 27.5, 25.0) == newton(3.0u"m", 30.0u"cm", 27.5u"cm", 25.0u"cm")

    # unit normalisation: results must stay coherent (m^3) when length and diameter
    # units differ, instead of leaking a composite unit like `cm m^2`
    @test unit(smalian(300.0u"cm", 30.0u"cm", 25.0u"cm")) == u"m^3"
    @test ustrip(uconvert(u"m^3", smalian(300.0u"cm", 30.0u"cm", 25.0u"cm"))) ≈ 0.17965982987716633 atol = 1e-6
    @test unit(huber(3000.0u"mm", 27.5u"cm")) == u"m^3"
    @test ustrip(uconvert(u"m^3", huber(3000.0u"mm", 27.5u"cm"))) ≈ 0.1781872083207961 atol = 1e-6
    @test unit(newton(10.0u"ft", 30.0u"cm", 27.5u"cm", 25.0u"cm")) == u"m^3"

    # imperial diameters keep results in cubic feet
    @test ustrip(uconvert(u"ft^3", smalian(10.0u"ft", 12.0u"inch", 10.0u"inch"))) ≈ 6.654067773228381 atol = 1e-6
    @test ustrip(uconvert(u"ft^3", huber(10.0u"ft", 11.0u"inch"))) ≈ 6.599526234103559 atol = 1e-6
    @test ustrip(uconvert(u"ft^3", newton(10.0u"ft", 12.0u"inch", 11.0u"inch", 10.0u"inch"))) ≈ 6.617706747145165 atol = 1e-6

    # a zero diameter at the tip of the bole is a legitimate profile point
    h_vec = [0.1, 1.3, 3.3, 5.3]u"m"
    d_vec = [30.0, 25.0, 18.0, 10.0]u"cm"
    @test ustrip(uconvert(u"m^3", smalian(h_vec, d_vec))) ≈ 0.17969909978533616 atol = 1e-6
    @test smalian([0.1, 1.3, 3.3, 5.3], [30.0, 25.0, 18.0, 10.0]) == smalian(h_vec, d_vec)

    L = 25.0u"m"
    d_closing = [30.0, 28.0, 25.0, 22.0, 19.0, 16.0, 13.0, 10.0, 7.0, 4.0, 0.0]u"cm"
    @test ustrip(uconvert(u"m^3", smalian(L, d_closing))) ≈ 0.6467753875577987 atol = 1e-6
    @test ustrip(uconvert(u"m^3", hohenadl(L, d_closing))) ≈ 0.6683024372181924 atol = 1e-6
    @test smalian(25.0, ustrip.(d_closing)) == smalian(L, d_closing)
    @test hohenadl(25.0, ustrip.(d_closing)) == hohenadl(L, d_closing)
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

    # unitless convenience: volume in m^3, heights in m, diameters in cm
    @test artificialformfactor(0.3378, 18.5, 22.7) == artificialformfactor(vt, ht, dbh)
    @test naturalformfactor(0.3378, 18.5, ustrip.(h_values), ustrip.(d_values)) == naturalformfactor(vt, ht, h_values, d_values)
    @test quotientform(18.5, 22.7, ustrip.(h_values), ustrip.(d_values)) == quotientform(ht, dbh, h_values, d_values)

    result_df = cubage(h_values, d_values; dlimit=2.8u"cm")
    @test isa(result_df, DataFrame)
    @test ustrip(uconvert(u"m^3", result_df.vt[1])) > 0

    # unitless call must match the unitful one exactly, including scalar kwargs
    result_df_bare = cubage(ustrip.(h_values), ustrip.(d_values); dlimit=2.8)
    @test all(ustrip.(collect(result_df[1, :])) .≈ ustrip.(collect(result_df_bare[1, :])))

    # kwargs may mix plain numbers with a unitful data vector
    result_df_mixed = cubage(h_values, d_values; dlimit=2.8)
    @test all(ustrip.(collect(result_df[1, :])) .≈ ustrip.(collect(result_df_mixed[1, :])))
  end

  @testset "Cubage over/under bark" begin
    h_values = [0.3, 1.3, 3.3, 5.3]u"m"
    d_values = [9.0, 7.0, 5.8, 5.1]u"cm"
    e_values = [1.2, 0.8, 0.6, 0.4]u"cm"

    result_df = cubage(h_values, d_values, e_values; ht=7.0u"m")
    @test isa(result_df, DataFrame)
    @test ustrip(uconvert(u"m^3", result_df.vtob[1])) > ustrip(uconvert(u"m^3", result_df.vtub[1]))

    # (h, d, e) has the same shape as the multi-tree (id, h, d), so it deliberately has
    # no unitless convenience (see the cubage docstring); h, d, e must carry units.
    @test_throws MethodError cubage(ustrip.(h_values), ustrip.(d_values), ustrip.(e_values); ht=7.0)
  end

  @testset "Multi-tree cubage" begin
    ids = [1, 1, 1, 1, 2, 2, 2, 2]
    h_values = [0.3, 1.3, 3.3, 5.3, 0.3, 1.3, 3.3, 5.3]u"m"
    d_values = [9.0, 7.0, 5.8, 5.1, 10.0, 8.0, 6.5, 5.5]u"cm"

    result_df = cubage(ids, h_values, d_values; ht=[7.0, 7.0, 7.0, 7.0, 8.0, 8.0, 8.0, 8.0])
    @test isa(result_df, DataFrame)
    @test result_df.id == [1, 2]
    @test nrow(result_df) == 2

    e_values = [1.2, 0.8, 0.6, 0.4, 1.2, 0.8, 0.6, 0.4]u"cm"
    result_df_bark = cubage(ids, h_values, d_values, e_values; ht=7.0)
    @test isa(result_df_bark, DataFrame)
    @test result_df_bark.id == [1, 2]

    # numeric tree ids are the common case; h and d must carry units here, since a
    # fully unitless (id, h, d) call is indistinguishable from the single-tree
    # over/under-bark cubage(h, d, e) triple. Neither has an automatic Real fallback
    # (see above), so this must raise a clean MethodError rather than silently
    # dispatching to the wrong function (id treated as heights).
    @test_throws MethodError cubage(ids, ustrip.(h_values), ustrip.(d_values); ht=7.0)
  end
end
