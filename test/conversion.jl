import Base: isapprox
import StaticArrays: FieldVector, SVector
@testset "Co-ordinate system conversion" begin
    Base.isapprox(a::T, b::T; kwargs...) where {T<:FieldVector} =
        all(isapprox(SVector(a), SVector(b); atol=sqrt(eps()), kwargs...))
    @testset "Fixed conversions" begin
        # ENU <-> XYZ
        @testset "bearing 0°" begin
            bearing = 0.
            let XYZfromENU = XYZfromENU(bearing),
                ENUfromXYZ = ENUfromXYZ(bearing)

                enus = [ENU(1., 0., 0.),ENU(1.,  1., 1.),ENU(1.,  1., 1.),ENU(1.,  1., 1.)]
                xyzs = [XYZ(0., -1, 0.),XYZ(1., -1., 1.),XYZ(1., -1., 1.),XYZ(1., -1., 1.)]
                for (enu, xyz) in zip(enus, xyzs)
                    @test XYZfromENU(enu) ≈ xyz
                    @test ENUfromXYZ(xyz) ≈ enu
                    @test XYZfromENU(ENUfromXYZ(xyz)) ≈ xyz
                    @test ENUfromXYZ(XYZfromENU(enu)) ≈ enu
                end
            end
        end

        # ENU <-> XYZ
        @testset "bearing -90°" begin
            bearing = -τ/4
            let XYZfromENU = XYZfromENU(bearing),
                ENUfromXYZ = ENUfromXYZ(bearing)

                enus = [ENU( 1., 0., 0.),ENU( 1.,  1., 1.)]
                xyzs = [XYZ(-1., 0., 0.),XYZ(-1., -1., 1.)]
                for (enu, xyz) in zip(enus, xyzs)
                    @test XYZfromENU(enu) ≈ xyz
                    @test ENUfromXYZ(xyz) ≈ enu
                    @test XYZfromENU(ENUfromXYZ(xyz)) ≈ xyz
                    @test ENUfromXYZ(XYZfromENU(enu)) ≈ enu
                end
            end

        end

        @testset "distance" begin
            @test euclidean_distance(XYZ(1., 1., 1.), XYZ(-1., -1., 1)) ≈ √(2)*2
        end
    end
end
