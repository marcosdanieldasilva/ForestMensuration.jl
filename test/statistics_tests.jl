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
end


