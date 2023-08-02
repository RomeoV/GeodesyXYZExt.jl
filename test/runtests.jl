using Geodesy
import Geodesy: ENU
using GeodesyXYZExt
using Test
import Tau: τ
import Unitful

################################################
### Helpers for testing approximate equality ###
################################################

@testset "GeodesyXYZExt" begin
    include("conversion.jl")
    include("unitful.jl")
end # @testset "GeodesyXYZExt"
