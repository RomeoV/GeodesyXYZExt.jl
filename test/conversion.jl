import Base: isapprox
import StaticArrays: FieldVector, SVector
@testset "Co-ordinate system conversion" begin
    Base.isapprox(a::T, b::T; kwargs...) where {T<:FieldVector} =
        all(isapprox(SVector(a), SVector(b); atol=sqrt(eps()), kwargs...))
    @testset "Fixed conversions" begin
        # ENU <-> XYZ
        @testset "bearing 0°" begin
            bearing = 0.
            let XYZfromENU_ = XYZfromENU(bearing),
                ENUfromXYZ_ = ENUfromXYZ(bearing)

                enus = [ENU(1., 0., 0.),ENU(1.,  1., 1.),ENU(1.,  1., 1.),ENU(1.,  1., 1.)]
                xyzs = [XYZ(0., -1, 0.),XYZ(1., -1., 1.),XYZ(1., -1., 1.),XYZ(1., -1., 1.)]
                for (enu, xyz) in zip(enus, xyzs)
                    @test XYZfromENU_(enu) ≈ xyz
                    @test ENUfromXYZ_(xyz) ≈ enu
                    @test XYZfromENU_(ENUfromXYZ_(xyz)) ≈ xyz
                    @test ENUfromXYZ_(XYZfromENU_(enu)) ≈ enu
                end
            end
        end

        # ENU <-> XYZ
        @testset "bearing -90°" begin
            bearing = -τ/4
            let XYZfromENU_ = XYZfromENU(bearing),
                ENUfromXYZ_ = ENUfromXYZ(bearing)

                enus = [ENU( 1., 0., 0.),ENU( 1.,  1., 1.)]
                xyzs = [XYZ(-1., 0., 0.),XYZ(-1., -1., 1.)]
                for (enu, xyz) in zip(enus, xyzs)
                    @test XYZfromENU_(enu) ≈ xyz
                    @test ENUfromXYZ_(xyz) ≈ enu
                    @test XYZfromENU_(ENUfromXYZ_(xyz)) ≈ xyz
                    @test ENUfromXYZ_(XYZfromENU_(enu)) ≈ enu
                end
            end

        end

        @testset "distance" begin
            @test euclidean_distance(XYZ(1., 1., 1.), XYZ(-1., -1., 1)) ≈ √(2)*2
        end
    end

    @testset "Different coordinate systems" begin
        @test XYZ(1., 2., 3.) + XYZ(1., 2., 3.) == XYZ(2., 4., 6.)
        @test ENU(1., 2., 3.) + ENU(1., 2., 3.) == ENU(2., 4., 6.)
        @test ECEF(1., 2., 3.) + ECEF(1., 2., 3.) == ECEF(2., 4., 6.)
        @test @fastmath XYZ(1., 2., 3.) + XYZ(1., 2., 3.) == XYZ(2., 4., 6.)
        @test @fastmath ENU(1., 2., 3.) + ENU(1., 2., 3.) == ENU(2., 4., 6.)
        @test @fastmath ECEF(1., 2., 3.) + ECEF(1., 2., 3.) == ECEF(2., 4., 6.)

        @test_throws CoordinateSystemError XYZ(1., 2., 3.) + ENU(1., 2., 3.)
        @test_throws CoordinateSystemError ENU(1., 2., 3.) + XYZ(1., 2., 3.)
        @test_throws CoordinateSystemError XYZ(1., 2., 3.) + ECEF(1., 2., 3.)
        @test_throws CoordinateSystemError @fastmath XYZ(1., 2., 3.) + ENU(1., 2., 3.)
        @test_throws CoordinateSystemError @fastmath ENU(1., 2., 3.) + XYZ(1., 2., 3.)
        @test_throws CoordinateSystemError @fastmath XYZ(1., 2., 3.) + ECEF(1., 2., 3.)
    end


    @testset "Readme example" begin
        DATUM=wgs84;
        a, b, c = 27.1127, 109.3497, 5.0
        origin_lla::LLA = LLA(a, b, c)
        bearing::Float64 = deg2rad(45)  # measured clockwise from north

        let XYZfromLLA = XYZfromLLA(origin_lla, bearing, DATUM),
            LLAfromXYZ = LLAfromXYZ(origin_lla, bearing, DATUM)

            origin::XYZ = XYZfromLLA(LLA(a, b, c))
            corners = XYZ[
                origin + XYZ(0., 10., 0),
                origin + XYZ(0., -10., 0)
            ]

            LLAfromXYZ.(corners)
            @test LLAfromXYZ(origin) ≈ origin_lla
        end
    end
end
