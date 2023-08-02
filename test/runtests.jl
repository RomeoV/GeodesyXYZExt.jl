using Geodesy
import Geodesy: ENU
using GeodesyXYZExt
using Test
import Tau: τ
import Unitful
import Base: isapprox
import StaticArrays: FieldVector, SVector

Base.isapprox(a::T, b::T; kwargs...) where {T<:FieldVector} =
    all(isapprox(SVector(a), SVector(b); atol=sqrt(eps()), kwargs...))

@testset "GeodesyXYZExt" begin
    include("conversion.jl")
    include("unitful.jl")
    include("readme_blocks.jl")
end # @testset "GeodesyXYZExt"
