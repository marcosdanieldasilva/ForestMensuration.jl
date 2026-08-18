@testset "horizontalpointsampling" begin
  data = DataFrame(
    point=[1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5],
    diameter=[25.0, 30.0, 20.0, 28.0, 22.0, 35.0, 30.0, 25.0, 20.0, 22.0, 30.0, 28.0, 26.0],
    volume=[0.35, 0.55, 0.22, 0.48, 0.28, 0.85, 0.55, 0.35, 0.22, 0.28, 0.55, 0.48, 0.40],
  )

  report = sampling(HorizontalPointSampling, :point, :diameter, :volume, 2.0, 6, 0.1, 10, data)
  @test report isa SamplingReport
  rt = resultTable(report)
  pt = pointTable(report)

  @test nrow(pt) == 5   # point 6 (zero trees) never appears as a row
  @test rt.n[1] == 6    # but is still counted in the sample size
  @test rt.N[1] == 100
  @test eltype(rt.vha) <: Unitful.Quantity

  @testset "Gha identity: EFi*gi sums to ntrees*baf exactly" begin
    baf = 2.0
    for row in eachrow(pt)
      @test ustrip(row.Gha) ≈ row.ntrees * baf atol = 1e-9
    end
  end

  @testset "units" begin
    reportU = sampling(HorizontalPointSampling, :point, :diameter, :volume, 2.0, 6, 0.1u"ha", 10u"ha", data)
    @test resultTable(reportU) == rt
  end

  @testset "zero-count points are not silently dropped" begin
    dataNoZero = DataFrame(
      point=[1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 5, 5, 5],
      diameter=data.diameter, volume=data.volume,
    )
    reportDropped = sampling(HorizontalPointSampling, :point, :diameter, :volume, 2.0, 5, 0.1, 10, dataNoZero)
    @test ustrip(resultTable(reportDropped).vha[1]) > ustrip(rt.vha[1])
  end

  @test_throws ArgumentError sampling(HorizontalPointSampling, :point, :diameter, :volume, 2.0, 4, 0.1, 10, data)
  @test_throws ArgumentError sampling(HorizontalPointSampling, :point, :diameter, :volume, -2.0, 6, 0.1, 10, data)
end
