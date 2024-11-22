@testset "Regression Function Tests" begin

  # Test 1: Insufficient Data Points (n < k + 2)
  @testset "Insufficient Data Points" begin
    # Create data with insufficient observations
    data_insufficient = DataFrame(y=[1.0], x=[2.0])
    try
      regression(:y, :x, data_insufficient)
      @test false  # Should not reach here
    catch e
      @test isa(e, ErrorException)
      @test occursin("There are not enough data points to perform regression", e.msg)
    end
  end

  # Test 2: Coercion Failure (Non-numeric data)
  @testset "Coercion Failure" begin
    data_non_numeric = DataFrame(y=["a", "b", "c", "d", "e", "f"], x=["a", "b", "c", "d", "e", "f"])
    try
      regression(:y, :x, data_non_numeric)
      @test false  # Should not reach here
    catch e
      @test isa(e, ErrorException)
      @test occursin("Unable to coerce variables", e.msg)
    end
  end

  # Test 3: Successful Regression
  @testset "Successful Regression" begin
    data_valid = DataFrame(y=[1.0, 2.0, 3.0, 4.0], x=[2.0, 3.0, 4.0, 5.0])
    models = regression(:y, :x, data_valid)
    @test length(models) > 0
    @test models[1] isa TableRegressionModel
  end


  @testset "Cylinder Volume Calculation Tests" begin
    # Test 1: Standard case
    h_standard = 18.5
    d_standard = 30.0
    expected_volume = 1.3076879420567515
    @test cylinder_volume(h_standard, d_standard) ≈ expected_volume atol = 1e-6

    # Test 2: Edge case with zero height (should throw an error)
    h_zero = 0.0
    @test_throws DomainError cylinder_volume(h_zero, d_standard)

    # Test 3: Edge case with zero diameter (should throw an error)
    d_zero = 0.0
    @test_throws DomainError cylinder_volume(h_standard, d_zero)

    # Test 4: Negative height (should throw an error)
    @test_throws DomainError cylinder_volume(-h_standard, d_standard)

    # Test 5: Negative diameter (should throw an error)
    @test_throws DomainError cylinder_volume(h_standard, -d_standard)
  end

end