@testset "dendrometry_and_inventory.jl" begin
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
    end

    @testset "DataFrame dmetrics Tests" begin
      dfDiameters = dmetrics(diameters) .|> ustrip
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
    end

    @testset "DataFrame hmetrics Tests" begin
      dfHeights = hmetrics(diameters, heights) .|> ustrip
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
    end

  end

  @testset "Frequency Tables Function Tests" begin
    # Test Data
    diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]
    species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"]
    data = DataFrame(species=species, diameters=diameters)

    # Test 1: Standard Case
    result_df = frequency_table(diameters, 2)
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
    @test_throws DomainError frequency_table(diameters, -2)

    # Test 3: Invalid Class Width (Zero)
    @test_throws DomainError frequency_table(diameters, 0)

    # Test 4: Verify Cumulative Frequencies
    @test result_df.Fi[end] == sum(result_df.fi)
    @test result_df.Fri[end] == 100.0

    # Test 5: Check if total frequency equals number of data points
    total_frequency = sum(result_df.fi)
    @test total_frequency == length(diameters)
  end

  @testset "Diametric Table Function Tests" begin
    # Test Data
    diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]
    species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"]
    data = DataFrame(species=species, diameters=diameters)

    # Test 1: Standard Case
    result_df = diametric_table(diameters, 2, plot_area=0.05)
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

    for col in names(expected_df)
      @test all(isapprox.(result_df[!, col], expected_df[!, col], atol=1e-5))
    end

    # Test 2: Invalid Class Width (Negative)
    @test_throws DomainError diametric_table(diameters, -2)

    # Test 3: Invalid Class Width (Zero)
    @test_throws DomainError diametric_table(diameters, 0)

    # Test 4: Invalid Plot Area (Negative)
    @test_throws DomainError diametric_table(diameters, 2, plot_area=-0.05)

    # Test 5: Diameters with Negative Values
    @test_throws DomainError diametric_table(-diameters, 2)

    # Test 6: Diameters with Zero Values
    @test_throws DomainError diametric_table([0, 10.5, 12.0], 2)

    # Test 7: Check if cumulative ng equals sum of ng
    @test result_df.∑ng[end] ≈ sum(result_df.ng) atol = 1e-5

    # Test 8: Check if cumulative ng_ha equals sum of ng_ha
    @test result_df.∑ng_ha[end] ≈ sum(result_df.ng_ha) atol = 1e-5
  end

end


