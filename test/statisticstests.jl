@testset "dendrometrics and distributiontables" begin
  @testset "Dendrometric Averages Function Tests" begin

    @testset "Basal Area Tests" begin
      dStandard = 30.0
      expectedBasalArea = 0.07068583470577035
      @test basalarea(dStandard) |> ustrip ≈ expectedBasalArea atol = 1e-6
      dZero = 0.0
      @test_throws DomainError basalarea(dZero)
      @test_throws DomainError basalarea(-dStandard)
    end

    diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]
    heights = [10.2, 11.5, 12.3, 14.1, 14.9, 16.5, 17.2, 18.0, 19.6, 21.2]
    plotArea = 0.05

    @testset "Diameter Metrics Tests" begin
      @test dm(diameters) |> ustrip ≈ 17.25 atol = 1e-4
      @test dg(diameters) |> ustrip ≈ 17.779904 atol = 1e-6
      @test dw(diameters) |> ustrip ≈ 18.6 atol = 1e-4
      @test dz(diameters) |> ustrip ≈ 17.266296 atol = 1e-6
      @test dd(diameters, plotArea) |> ustrip ≈ 21.0 atol = 1e-4
      hohenadl = dh(diameters)
      @test hohenadl.dl |> ustrip ≈ 12.708524 atol = 1e-6
      @test hohenadl.du |> ustrip ≈ 21.791475 atol = 1e-6

      # a plot too small for the dominant-tree standard returns NaN, but the result
      # must still carry the diameter unit so DataFrame columns stay type-consistent
      dNaN = dd(diameters, 0.0)
      @test isnan(dNaN)
      @test unit(dNaN) == u"cm"
      @test unit(dd(diameters * u"inch", 0.0u"ac")) == u"inch"
    end

    @testset "DataFrame dmetrics Tests" begin
      dfDiametersRaw = dmetrics(diameters)
      @test isnan(dfDiametersRaw.dd[1])
      @test unit(dfDiametersRaw.dd[1]) == u"cm"

      dfDiameters = dfDiametersRaw .|> ustrip
      @test dfDiameters.dl[1] ≈ 12.708524 atol = 1e-6
      @test dfDiameters.dm[1] ≈ 17.25 atol = 1e-4
      @test dfDiameters.dg[1] ≈ 17.779904 atol = 1e-6
      @test dfDiameters.dw[1] ≈ 18.6 atol = 1e-4
      @test dfDiameters.dz[1] ≈ 17.266296 atol = 1e-6
      @test isnan(dfDiameters.dd[1])
      @test dfDiameters.du[1] ≈ 21.791475 atol = 1e-6
      @test dfDiameters.dv[1] ≈ 26.3274 atol = 1e-4
      dfDiametersArea = dmetrics(diameters, plotArea)
      @test dfDiametersArea.dd[1] |> ustrip ≈ 21.0 atol = 1e-4
    end

    @testset "Height Metrics Tests" begin
      @test hm(heights) |> ustrip ≈ 15.55 atol = 1e-4
      @test hd(diameters, heights, plotArea) |> ustrip ≈ 18.5 atol = 1e-4
      @test hg(diameters, heights) |> ustrip ≈ 17.148042 atol = 1e-6

      hNaN = hd(diameters, heights, 0.0)
      @test isnan(hNaN)
      @test unit(hNaN) == u"m"
    end

    @testset "DataFrame hmetrics Tests" begin
      dfHeightsRaw = hmetrics(diameters, heights)
      @test isnan(dfHeightsRaw.hd[1])
      @test unit(dfHeightsRaw.hd[1]) == u"m"

      dfHeights = dfHeightsRaw .|> ustrip
      @test dfHeights.hl[1] ≈ 11.9589 atol = 1e-4
      @test dfHeights.hm[1] ≈ 15.55 atol = 1e-4
      @test isnan(dfHeights.hd[1])
      @test dfHeights.hg[1] ≈ 17.1480 atol = 1e-4
      @test dfHeights.hu[1] ≈ 19.1411 atol = 1e-4
      @test dfHeights.hv[1] ≈ 23.094 atol = 1e-3
      dfHeightsArea = hmetrics(diameters, heights, plotArea)
      @test dfHeightsArea.hd[1] |> ustrip ≈ 18.5 atol = 1e-4
    end

    @testset "DataFrame standmetrics Tests" begin
      dfStand = standmetrics(diameters, heights, plotArea) .|> ustrip
      @test dfStand.n[1] == 10
      @test dfStand.g[1] ≈ 0.248284 atol = 1e-6
      @test dfStand.EF[1] ≈ 20.0 atol = 1e-4
      @test dfStand.N[1] ≈ 200.0 atol = 1e-4
      @test dfStand.G[1] ≈ 4.96568 atol = 1e-5
      @test dfStand.dl[1] ≈ 12.7085 atol = 1e-4
      @test dfStand.dm[1] ≈ 17.25 atol = 1e-4
      @test dfStand.dg[1] ≈ 17.7799 atol = 1e-4
      @test dfStand.dw[1] ≈ 18.6 atol = 1e-4
      @test dfStand.dz[1] ≈ 17.2663 atol = 1e-4
      @test dfStand.dd[1] ≈ 21.0 atol = 1e-4
      @test dfStand.du[1] ≈ 21.7915 atol = 1e-4
      @test dfStand.dv[1] ≈ 26.3274 atol = 1e-4
      @test dfStand.hl[1] ≈ 11.9589 atol = 1e-4
      @test dfStand.hm[1] ≈ 15.55 atol = 1e-4
      @test dfStand.hd[1] ≈ 18.5 atol = 1e-4
      @test dfStand.hg[1] ≈ 17.148 atol = 1e-3
      @test dfStand.hu[1] ≈ 19.1411 atol = 1e-4
      @test dfStand.hv[1] ≈ 23.094 atol = 1e-3

      # results follow the unit of the input diameters/heights, not a hardcoded one
      dfStandInch = standmetrics(diameters * u"inch", heights * u"ft", plotArea * u"ac")
      @test unit(dfStandInch.dm[1]) == u"inch"
      @test unit(dfStandInch.hm[1]) == u"ft"
      @test unit(dfStandInch.EF[1]) == u"ac^-1"
    end

  end

  @testset "Frequency Tables Function Tests" begin
    # Test Data
    diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]
    species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"]
    data = DataFrame(species=species, diameters=diameters)

    # Test 1: Standard Case
    result_df = frequencytable(diameters, 2)
    @test isa(result_df, DataFrame)

    # Expected results for the test data
    expected_df = DataFrame(
      LI=[10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0],
      Xi=[11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0],
      LS=[12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0],
      fi=[1, 1, 1, 2, 1, 2, 1, 1],
      Fi=[1, 2, 3, 5, 6, 8, 9, 10],
      fri=[10.0, 10.0, 10.0, 20.0, 10.0, 20.0, 10.0, 10.0],
      Fri=[10.0, 20.0, 30.0, 50.0, 60.0, 80.0, 90.0, 100.0]
    )

    # Compare the result with expected values
    @test size(result_df) == size(expected_df)

    for col in names(expected_df)
      @test result_df[!, col] == expected_df[!, col]
    end

    # Test 2: Invalid Class Width (Negative)
    @test_throws DomainError frequencytable(diameters, -2)

    # Test 3: Invalid Class Width (Zero)
    @test_throws DomainError frequencytable(diameters, 0)

    # Test 4: Verify Cumulative Frequencies
    @test result_df.Fi[end] == sum(result_df.fi)
    @test result_df.Fri[end] == 100.0

    # Test 5: Check if total frequency equals number of data points
    total_frequency = sum(result_df.fi)
    @test total_frequency == length(diameters)

    @testset "with units" begin
      # a unitful sample keeps its unit in the class limits, and a bare class width
      # is interpreted in that same unit
      result_df_u = frequencytable(diameters * u"cm", 2u"cm")
      @test ustrip.(result_df_u.LI) == expected_df.LI
      @test ustrip.(result_df_u.Xi) == expected_df.Xi
      @test ustrip.(result_df_u.LS) == expected_df.LS
      @test unit(result_df_u.LI[1]) == u"cm"

      result_df_u2 = frequencytable(diameters * u"cm", 2)
      @test result_df_u2 == result_df_u

      # frequencytable is not diameter-specific: any unitful variable works, e.g. heights
      result_df_h = frequencytable(diameters * u"m", 2u"m")
      @test unit(result_df_h.LI[1]) == u"m"
      @test ustrip.(result_df_h.LI) == expected_df.LI

      # auto class width also respects the sample's unit
      auto_u = frequencytable(diameters * u"cm")
      auto_plain = frequencytable(diameters)
      @test ustrip.(auto_u.LI) == auto_plain.LI
      @test unit(auto_u.LI[1]) == u"cm"
    end

    @testset "grouped" begin
      result_df_g = frequencytable(:species, :diameters, data)
      @test isa(result_df_g, DataFrame)

      dataU = DataFrame(species=species, diameters=diameters * u"cm")
      result_df_gu = frequencytable(:species, :diameters, dataU)
      @test ustrip.(result_df_gu.LI) == result_df_g.LI
      @test unit(result_df_gu.LI[1]) == u"cm"

      result_df_gu_hi = frequencytable(:species, :diameters, 2u"cm", dataU)
      @test unit(result_df_gu_hi.LI[1]) == u"cm"
    end
  end

  @testset "Diametric Table Function Tests" begin
    # Test Data
    diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]
    species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"]
    data = DataFrame(species=species, diameters=diameters)

    # Test 1: Standard Case
    result_df = diametrictable(diameters, 2, plot_area=0.05)
    @test isa(result_df, DataFrame)

    # Expected results for the test data
    expected_df = DataFrame(
      LI=[10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0],
      Xi=[11.0, 13.0, 15.0, 17.0, 19.0, 21.0, 23.0, 25.0],
      LS=[12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0],
      fi=[1, 1, 1, 2, 1, 2, 1, 1],
      Fi=[1, 2, 3, 5, 6, 8, 9, 10],
      fri=[10.0, 10.0, 10.0, 20.0, 10.0, 20.0, 10.0, 10.0],
      Fri=[10.0, 20.0, 30.0, 50.0, 60.0, 80.0, 90.0, 100.0],
      g=[0.00950332, 0.0132732, 0.0176715, 0.022698, 0.0283529, 0.0346361, 0.0415476, 0.0490874],
      ng=[0.00950332, 0.0132732, 0.0176715, 0.045396, 0.0283529, 0.0692721, 0.0415476, 0.0490874],
      ∑ng=[0.00950332, 0.0227765, 0.040448, 0.085844, 0.114197, 0.183469, 0.225017, 0.274104],
      fi_ha=[20.0, 20.0, 20.0, 40.0, 20.0, 40.0, 20.0, 20.0],
      Fi_ha=[20.0, 40.0, 60.0, 100.0, 120.0, 160.0, 180.0, 200.0],
      ng_ha=[0.190066, 0.265465, 0.353429, 0.90792, 0.567057, 1.38544, 0.830951, 0.981748],
      ∑ng_ha=[0.190066, 0.455531, 0.80896, 1.71688, 2.28394, 3.66938, 4.50033, 5.48208]
    )

    # Compare the result with expected values
    @test size(result_df) == size(expected_df)

    # like dm/dg/... elsewhere in the package, diametrictable always returns unitful
    # quantities even when called with plain numbers, so strip units before comparing
    for col in names(expected_df)
      @test all(isapprox.(ustrip.(result_df[!, col]), expected_df[!, col], atol=1e-5))
    end

    # Test 2: Invalid Class Width (Negative)
    @test_throws DomainError diametrictable(diameters, -2)

    # Test 3: Invalid Class Width (Zero)
    @test_throws DomainError diametrictable(diameters, 0)

    # Test 4: Invalid Plot Area (Negative)
    @test_throws DomainError diametrictable(diameters, 2, plot_area=-0.05)

    # Test 5: Diameters with Negative Values
    @test_throws DomainError diametrictable(-diameters, 2)

    # Test 6: Diameters with Zero Values
    @test_throws DomainError diametrictable([0, 10.5, 12.0], 2)

    # Test 7: Check if cumulative ng equals sum of ng
    @test ustrip(result_df.∑ng[end]) ≈ ustrip(sum(result_df.ng)) atol = 1e-5

    # Test 8: Check if cumulative ng_ha equals sum of ng_ha
    @test ustrip(result_df.∑ng_ha[end]) ≈ ustrip(sum(result_df.ng_ha)) atol = 1e-5

    @testset "with units" begin
      # diameters and class width as quantities, plot area as a quantity in m^2
      result_df_u = diametrictable(diameters * u"cm", 2u"cm", plot_area=500u"m^2")
      for col in names(expected_df)
        @test all(isapprox.(ustrip.(result_df_u[!, col]), expected_df[!, col], atol=1e-5))
      end
      @test unit(result_df_u.LI[1]) == u"cm"
      @test unit(result_df_u.g[1]) == u"m^2"
      @test unit(result_df_u.fi_ha[1]) == u"ha^-1"

      # mixing a unitful diameter vector with a bare class width and plot area
      result_df_mixed = diametrictable(diameters * u"cm", 2, plot_area=0.05)
      @test all(isapprox.(ustrip.(result_df_mixed[!, :g]), expected_df.g, atol=1e-5))

      # a plot exactly one reference area wide needs no expansion, so those columns
      # are omitted instead of being a redundant copy of fi/ng
      result_df_default = diametrictable(diameters, 2)
      @test "fi_ha" ∉ names(result_df_default)
      @test "Fi_ha" ∉ names(result_df_default)
      @test "ng_ha" ∉ names(result_df_default)
      @test "∑ng_ha" ∉ names(result_df_default)
      @test ncol(result_df_default) == ncol(result_df) - 4

      # imperial diameters expand per acre and report areas in ft^2
      imperial = diametrictable([4.0, 5.0, 6.0, 7.0]u"inch", 1u"inch", plot_area=0.1u"ac")
      @test unit(imperial.g[1]) == u"ft^2"
      @test unit(imperial.fi_ha[1]) == u"ac^-1"
    end

    @testset "grouped with units" begin
      dataU = DataFrame(species=species, diameters=diameters * u"cm")
      result_df_gu = diametrictable(:species, :diameters, dataU, plot_area=0.05u"ha")
      @test unit(result_df_gu.LI[1]) == u"cm"
      @test unit(result_df_gu.g[1]) == u"m^2"

      result_df_gu_hi = diametrictable(:species, :diameters, 3u"cm", dataU, plot_area=0.05u"ha")
      @test unit(result_df_gu_hi.LI[1]) == u"cm"
    end
  end

end
