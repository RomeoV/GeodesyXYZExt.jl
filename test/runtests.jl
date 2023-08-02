using Geodesy
import Geodesy: ENU
using GeodesyXYZExt
using Test
import Tau: τ

################################################
### Helpers for testing approximate equality ###
################################################

# Interesting that this isn't in Base...
Base.isapprox(a::T, b::T; kwargs...) where {T<:Tuple} = all(ntuple(i->isapprox(a[i],b[i]), length(a)); kwargs...)

@testset "Geodesy" begin
    include("conversion.jl")
end # @testset "GeodesyXYZExt"
