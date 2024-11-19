@testset "cubage.jl" begin
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
    @test_throws DomainError cylinder_volume(h_zero, d_standard)

    # Test 3: Edge case with zero diameter (should throw an error)
    d_zero = 0.0
    @test_throws DomainError cylinder_volume(h_standard, d_zero)

    # Test 4: Negative height (should throw an error)
    @test_throws DomainError cylinder_volume(-h_standard, d_standard)

    # Test 5: Negative diameter (should throw an error)
    @test_throws DomainError cylinder_volume(h_standard, -d_standard)
  end

  @testset "Cone Volume Calculation Tests" begin
    # Test 1: Standard case
    h_standard = 18.5
    d_standard = 30.0
    expected_volume = 0.4358959806855838
    @test cone_volume(h_standard, d_standard) ≈ expected_volume atol = 1e-6

    # Test 2: Edge case with zero height (should throw an error)
    h_zero = 0.0
    @test_throws DomainError cone_volume(h_zero, d_standard)

    # Test 3: Edge case with zero diameter (should throw an error)
    d_zero = 0.0
    @test_throws DomainError cone_volume(h_standard, d_zero)

    # Test 4: Negative height (should throw an error)
    @test_throws DomainError cone_volume(-h_standard, d_standard)

    # Test 5: Negative diameter (should throw an error)
    @test_throws DomainError cone_volume(h_standard, -d_standard)
  end

  @testset "Bark Factor Calculation Tests" begin
    # Test 1: Standard case with positive diameters and non-negative bark thicknesses
    d_values = [30.0, 22.5, 20.2, 15.4, 13.2, 10.9]
    e_values = [1.2, 1.1, 0.85, 0.66, 0.48, 0.0]  # Bark thicknesses
    expected_factor = 0.961764705882353
    @test bark_factor(d_values, e_values) ≈ expected_factor atol = 1e-6

    # Test 2: Error case with a negative diameter (should throw an error)
    @test_throws DomainError bark_factor(-d_values, e_values)

    # Test 3: Error case with a negative bark thickness (should throw an error)
    @test_throws DomainError bark_factor(d_values, -e_values)
  end

  @testset "Bole Volume Function Tests" begin
    d_values = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0]
    h_values = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8]

    # Test 1: Smalian method
    expected_volume_smalian = 0.021087744337680632
    @test bole_volume(Smalian, h_values, d_values) ≈ expected_volume_smalian atol = 1e-6

    # Test 2: Huber method
    expected_volume_huber = 0.020708986073382216
    @test bole_volume(Huber, h_values, d_values) ≈ expected_volume_huber atol = 1e-6

    # Test 3: Newton method
    expected_volume_newton = 0.015548265641391484
    @test bole_volume(Newton, h_values, d_values) ≈ expected_volume_newton atol = 1e-6

    # Test 4: Error when h and d have different lengths
    h_values_short = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3]
    @test_throws DimensionMismatch bole_volume(Smalian, h_values_short, d_values)

    # Test 5: Error when h contains negative values
    @test_throws DomainError bole_volume(Smalian, -h_values, d_values)

    # Test 6: Error when d contains negative values
    @test_throws DomainError bole_volume(Smalian, h_values, -d_values)
  end

  # Test set for Artificial Form Factor
  @testset "Artificial Form Factor Function Tests" begin
    vt = 0.3378
    ht = 18.5
    dbh = 22.7

    # Test 1: Standard case
    expected_aff = 0.451176344374475
    @test artificial_form_factor(vt, ht, dbh) ≈ expected_aff atol = 1e-6

    # Test 2: Error when vt <= 0
    @test_throws DomainError artificial_form_factor(0.0, ht, dbh)
    @test_throws DomainError artificial_form_factor(-1.0, ht, dbh)
  end

  # Test set for Natural Form Factor
  @testset "Natural Form Factor Function Tests" begin
    vt = 0.3378
    ht = 18.5
    d_values = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0]
    h_values = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8]

    # Test 1: Standard case
    expected_nff = 8.951469588617691
    @test natural_form_factor(vt, ht, h_values, d_values) ≈ expected_nff atol = 1e-6

    # Test 2: Error when vt <= 0
    @test_throws DomainError natural_form_factor(0.0, ht, h_values, d_values)
    @test_throws DomainError natural_form_factor(-1.0, ht, h_values, d_values)
  end

  @testset "Form Quotient Function Tests" begin
    vt = 0.3378
    ht = 18.5
    d_values = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0]
    h_values = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8]

    # Test 1: Standard case
    expected_nff = 8.951469588617691
    @test natural_form_factor(vt, ht, h_values, d_values) ≈ expected_nff atol = 1e-6

    # Test 2: Error when vt <= 0
    @test_throws DomainError natural_form_factor(0.0, ht, h_values, d_values)
    @test_throws DomainError natural_form_factor(-1.0, ht, h_values, d_values)
  end

  @testset "Cubage Function Tests" begin
    d_values = [9.0, 7.0, 5.8, 5.1, 3.8, 1.9, 0.0]
    h_values = [0.3, 1.3, 3.3, 5.3, 7.3, 9.3, 10.8]
    d_limit = 2.8

    # Test 1: Standard case
    result_df = cubage(Smalian, h_values, d_values, d_limit)
    @test isa(result_df, DataFrame)

    expected_values = Dict(
      :vt => 0.0228547,
      :v0 => 0.00190852,
      :vc => 0.0203784,
      :vr => 0.000425975,
      :vn => 0.000141764,
      :dbh => 7.0,
      :ht => 10.8,
      :hc => 8.35263,
      :aff => 0.549877,
      :nff => 0.486761,
      :qf => 0.719286
    )

    for (key, expected) in expected_values
      @test result_df[!, key][1] ≈ expected atol = 1e-4
    end

    # Test 2: Error when d_limit is negative
    @test_throws DomainError cubage(Smalian, h_values, d_values, -1.0)

    # Test 3: Error when dbh is not positive
    @test_throws DomainError cubage(Smalian, h_values, d_values, d_limit; dbh=0.0)
    @test_throws DomainError cubage(Smalian, h_values, d_values, d_limit; dbh=-1.0)
  end
end