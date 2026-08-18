@testset "clustersampling" begin
  data = DataFrame(
    cluster=repeat(1:6, inner=4),
    volume=[18.2, 19.1, 17.8, 18.9, 22.4, 23.1, 21.9, 22.8, 15.1, 16.0, 14.8, 15.6,
      27.3, 28.1, 26.9, 27.8, 19.8, 20.5, 19.1, 20.0, 24.5, 25.2, 23.9, 24.8],
  )

  report = clustersampling(:cluster, :volume, 0.02, 15, data)
  @test report isa SamplingReport
  rt = report.result_table

  @test ustrip(rt.vm[1]) ≈ 21.4 atol = 1e-6
  @test rt.M[1] == 4
  @test rt.n[1] == 6
  @test rt.N[1] == 188
  @test rt.icc[1] ≈ 0.9843 atol = 1e-3
  @test eltype(rt.vm) <: Unitful.Quantity

  @testset "M=1 matches simplecasualsampling" begin
    v = [381.7, 458.9, 468.2, 531.7, 474.1, 401.9, 469.1, 437.4, 435.3, 403.2, 397.1]
    data1 = DataFrame(cluster=1:11, volume=v)
    rc = clustersampling(:cluster, :volume, 0.05, 10, data1).result_table
    rs = simplecasualsampling(v, 0.05, 10)
    @test ustrip(rc.vm[1]) ≈ ustrip(rs.vm[1]) atol = 1e-6
    @test rc.nreq[1] == rs.nreq[1]
    @test rc.icc[1] ≈ 1.0 atol = 1e-9
  end

  @testset "units" begin
    reportU = clustersampling(:cluster, :volume, 0.02u"ha", 15u"ha", data)
    @test reportU.result_table == report.result_table
  end

  @test_throws ArgumentError clustersampling(:cluster, :volume, 0.02, 15,
    DataFrame(cluster=[1, 1, 2, 2, 2], volume=[1.0, 2.0, 3.0, 4.0, 5.0]))
end
