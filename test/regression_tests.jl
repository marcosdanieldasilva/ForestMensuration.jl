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


  @testset "Site Classifications Tests" begin
    data_age = DataFrame(
      age=repeat([36, 48, 60, 72, 84], outer=6),
      h=[13.6, 17.8, 21.5, 21.5, 21.8,
        14.3, 17.8, 21.0, 21.0, 21.4,
        14.0, 17.5, 21.2, 21.2, 21.4,
        13.4, 18.0, 20.8, 20.8, 23.2,
        13.2, 17.4, 20.3, 20.3, 22.0,
        13.2, 17.8, 21.3, 21.3, 22.5]
    )

    index_age = 60

    # Test 1: Standard case
    reg = regression(:h, :age, data_age) |> criteria_selection

    expected_site = [
      20.5, 20.3, 21.5, 20.4, 20.4, 22.2, 20.3, 21.0, 19.9, 20.0,
      21.4, 19.9, 21.2, 20.1, 20.0, 20.0, 20.5, 20.8, 19.8, 21.6,
      19.5, 19.7, 20.3, 19.3, 20.5, 19.5, 20.3, 21.3, 20.2, 21.0]

    @test site_classification(reg, data_age, index_age) ≈ expected_site atol = 1e-6
    @test site_classification(reg, index_age) ≈ expected_site atol = 1e-6

    # Test 2: Edge case with zero data age (should throw an error)
    index_age_zero = 0.0
    @test_throws DomainError site_classification(reg, data_age, index_age_zero)
    @test_throws DomainError site_classification(reg, index_age_zero)

    # Test 3: Edge case with negative data age (should throw an error)
    @test_throws DomainError site_classification(reg, data_age, -index_age)
    @test_throws DomainError site_classification(reg, -index_age)

    expected_height = [
      13.6, 17.8, 21.5, 21.5, 21.8, 14.3, 17.8, 21.0, 21.0, 21.4,
      14.0, 17.5, 21.2, 21.2, 21.4, 13.4, 18.0, 20.8, 20.9, 23.2,
      13.2, 17.4, 20.3, 20.3, 22.0, 13.2, 17.8, 21.3, 21.3, 22.5
    ]

    @test hdom_classification(reg, data_age, index_age, expected_site) ≈ expected_height atol = 1e-6

    hi = 1
    hi_zero = 0.0

    # Test 3: Edge case with negative data age (should throw an error)
    @test_throws DomainError site_table(reg, index_age, -hi)
    @test_throws DomainError site_table(reg, index_age, hi_zero)
  end

end