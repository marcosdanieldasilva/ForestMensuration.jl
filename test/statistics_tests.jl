@testset "dendrometry.jl" begin
  @testset "Basal Area Tests" begin
    # Test 1: Standard case
    d_standard = 30.0
    expected_basal_area = 0.07068583470577035
    @test basal_area(d_standard) ≈ expected_basal_area atol = 1e-6

    # Test 2: Edge case with zero diameter (should throw an error)
    d_zero = 0.0
    @test_throws DomainError basal_area(d_zero)

    # Test 3: Negative diameter (should throw an error)
    @test_throws DomainError basal_area(-d_standard)
  end

  @testset "Dendrometric Averages Function Tests" begin
    # Test Data
    diameters = [10.5, 12.0, 13.5, 15.0, 16.5, 18.0, 19.5, 21.0, 22.5, 24.0]
    species = ["Oak", "Oak", "Oak", "Oak", "Oak", "Pine", "Pine", "Pine", "Pine", "Pine"]
    data = DataFrame(species=species, diameters=diameters)

    # Test 1: Standard Case
    result_df = dendrometric_averages(diameters, area=0.05)
    @test isa(result_df, DataFrame)

    # Expected results
    expected_values = Dict(
      :d₋ => 12.7085,
      :d̅ => 17.25,
      :dg => 17.7799,
      :dw => 18.6,
      :dz => 17.2663,
      :d₁₀₀ => 21.0,
      :d₊ => 21.7915
    )

    # Compare the calculated values with expected values
    for (key, expected) in expected_values
      @test result_df[!, key][1] ≈ expected atol = 1e-4
    end

    # Test 2: Zero or Negative Area Value
    @test_throws DomainError dendrometric_averages(diameters, area=-1.0)
    @test_throws DomainError dendrometric_averages(diameters, area=0.0)

    # Test 3: Zero or Negative Diameters

    @test_throws DomainError dendrometric_averages([0, 10.5, 12.0])
    @test_throws DomainError dendrometric_averages(-diameters)

    # Test 4: Grouped Dendrometric Averages
    result_df = dendrometric_averages(:species, :diameters, data; area=0.02)
    @test isa(result_df, DataFrame)

    # Expected results for the test data
    expected_results = DataFrame(
      species=["Oak", "Pine"],
      d₋=[11.1283, 18.6283],
      d̅=[13.5, 21.0],
      dg=[13.6657, 21.1069],
      dw=[14.1, 21.6],
      dz=[13.5, 21.0],
      d₁₀₀=[15.75, 23.25],
      d₊=[15.8717, 23.3717]
    )

    for col in names(expected_results)[2:end]  # Skip 'species' column
      @test result_df[!, col] ≈ expected_results[!, col] atol = 1e-4
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


