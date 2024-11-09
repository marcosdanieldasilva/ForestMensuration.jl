using ForestMensuration
using Test

# Define shared values for tests
@testset "Diameter Interpolation Tests" begin
  h_values = [1.0, 1.3, 2.0, 3.0, 4.0, 5.0]
  d_values = [30.0, 22.5, 20.2, 15.4, 13.2, 10.9]

  # Test 1: Standard interpolation within the range
  h0 = 2.5  # Height between 2.0 and 3.0
  expected_diameter = 17.8
  @test ForestMensuration._diameter_interpolation(h0, h_values, d_values) ≈ expected_diameter atol = 1e-4

  # Test 2: Edge case where h0 exactly matches an element in h
  h0_exact = 3.0
  expected_diameter_exact = 15.4
  @test ForestMensuration._diameter_interpolation(h0_exact, h_values, d_values) ≈ expected_diameter_exact atol = 1e-4

  # Test 3: Edge case where h0 is exactly at the lower bound of h
  h0_lower = 1.0
  expected_diameter_lower = 30.0
  @test ForestMensuration._diameter_interpolation(h0_lower, h_values, d_values) ≈ expected_diameter_lower atol = 1e-4

  # Test 4: Edge case where h0 is exactly at the upper bound of h
  h0_upper = 5.0
  expected_diameter_upper = 10.9
  @test ForestMensuration._diameter_interpolation(h0_upper, h_values, d_values) ≈ expected_diameter_upper atol = 1e-4

  # Test 5: Error case where h0 is outside the range of h (below minimum)
  h0_below = 0.5
  @test_throws ErrorException ForestMensuration._diameter_interpolation(h0_below, h_values, d_values)

  # Test 6: Error case where h0 is outside the range of h (above maximum)
  h0_above = 6.0
  @test_throws ErrorException ForestMensuration._diameter_interpolation(h0_above, h_values, d_values)
end

@testset "Height Interpolation Tests" begin
  h_values = [1.0, 1.3, 2.0, 3.0, 4.0, 5.0]
  d_values = [30.0, 22.5, 20.2, 15.4, 13.2, 10.9]

  # Test 1: Standard interpolation within the range
  d_limit = 17.8  # Diameter between 20.2 and 15.4
  expected_height = 2.5
  expected_index = 4
  interpolated_height, index = ForestMensuration._height_interpolation(d_limit, h_values, d_values)
  @test interpolated_height ≈ expected_height atol = 1e-4
  @test index == expected_index

  # Test 2: Edge case where d_limit exactly matches an element in d
  d_limit_exact = 15.4
  expected_height_exact = 3.0
  expected_index_exact = 4
  interpolated_height, index = ForestMensuration._height_interpolation(d_limit_exact, h_values, d_values)
  @test interpolated_height ≈ expected_height_exact atol = 1e-4
  @test index == expected_index_exact

  # Test 3: Edge case where d_limit is exactly at the upper bound of d
  d_limit_upper = 30.0
  expected_height_upper = 1.0
  expected_index_upper = 2
  interpolated_height, index = ForestMensuration._height_interpolation(d_limit_upper, h_values, d_values)
  @test interpolated_height ≈ expected_height_upper atol = 1e-4
  @test index == expected_index_upper

  # Test 4: Edge case where d_limit is exactly at the lower bound of d
  d_limit_lower = 10.9
  expected_height_lower = 5.0
  expected_index_lower = 6
  interpolated_height, index = ForestMensuration._height_interpolation(d_limit_lower, h_values, d_values)
  @test interpolated_height ≈ expected_height_lower atol = 1e-4
  @test index == expected_index_lower

  # Test 5: Error case where d_limit is outside the range of d (above maximum)
  d_limit_above = 35.0
  @test_throws ErrorException ForestMensuration._height_interpolation(d_limit_above, h_values, d_values)

  # Test 6: Error case where d_limit is outside the range of d (below minimum)
  d_limit_below = 5.0
  @test_throws ErrorException ForestMensuration._height_interpolation(d_limit_below, h_values, d_values)
end

@testset "Cylinder Volume Calculation Tests" begin
  # Test 1: Standard case
  h_standard = 18.5
  d_standard = 30.0
  expected_volume = 1.3076879420567515
  @test cylinder_volume(h_standard, d_standard) ≈ expected_volume atol = 1e-6

  # Test 2: Edge case with zero height (should throw an error)
  h_zero = 0.0
  @test_throws ArgumentError cylinder_volume(h_zero, d_standard)

  # Test 3: Edge case with zero diameter (should throw an error)
  d_zero = 0.0
  @test_throws ArgumentError cylinder_volume(h_standard, d_zero)

  # Test 4: Negative height (should throw an error)
  @test_throws ArgumentError cylinder_volume(-h_standard, d_standard)

  # Test 5: Negative diameter (should throw an error)
  @test_throws ArgumentError cylinder_volume(h_standard, -d_standard)
end

@testset "Cone Volume Calculation Tests" begin
  # Test 1: Standard case
  h_standard = 18.5
  d_standard = 30.0
  expected_volume = 0.4358959806855838
  @test cone_volume(h_standard, d_standard) ≈ expected_volume atol = 1e-6

  # Test 2: Edge case with zero height (should throw an error)
  h_zero = 0.0
  @test_throws ArgumentError cone_volume(h_zero, d_standard)

  # Test 3: Edge case with zero diameter (should throw an error)
  d_zero = 0.0
  @test_throws ArgumentError cone_volume(h_standard, d_zero)

  # Test 4: Negative height (should throw an error)
  @test_throws ArgumentError cone_volume(-h_standard, d_standard)

  # Test 5: Negative diameter (should throw an error)
  @test_throws ArgumentError cone_volume(h_standard, -d_standard)
end

@testset "Bark Factor Calculation Tests" begin
  # Test 1: Standard case with positive diameters and non-negative bark thicknesses
  d_values = [30.0, 22.5, 20.2, 15.4, 13.2, 10.9]
  e_values = [1.2, 1.1, 0.85, 0.66, 0.48, 0.0]  # Bark thicknesses
  expected_factor = 0.961764705882353
  @test bark_factor(d_values, e_values) ≈ expected_factor atol = 1e-6

  # Test 2: Error case with a negative diameter (should throw an error)
  @test_throws ArgumentError bark_factor(-d_values, e_values)

  # Test 3: Error case with a negative bark thickness (should throw an error)
  @test_throws ArgumentError bark_factor(d_values, -e_values)
end
