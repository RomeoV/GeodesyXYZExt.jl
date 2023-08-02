@testset "Readme.md blocks" begin
    @testset "block 1" begin
        using GeodesyXYZExt; using Geodesy
        DATUM=wgs84;
        a, b, c = 27.1127, 109.3497, 5.0
        origin_lla::LLA = LLA(a, b, c)
        bearing::Float64 = deg2rad(45)  # measured clockwise from north

        # let's bind the transformations for easier handling
        let XYZfromLLA_ = XYZfromLLA(origin_lla, bearing, DATUM),
            LLAfromXYZ_ = LLAfromXYZ(origin_lla, bearing, DATUM)

            origin::XYZ = XYZfromLLA_(LLA(a, b, c))
            corners = XYZ[
                origin + XYZ(0., 10., 0),
                origin + XYZ(0., -10., 0)
            ]

            LLAfromXYZ_.(corners)
            @test LLAfromXYZ_(origin) ≈ origin_lla
        end
    end
    @testset "block 2" begin
        using GeodesyXYZExt
        GeodesyXYZExt.fixdatum!(wgs84)
        GeodesyXYZExt.fixorigin!(LLA(27.1127, 109.3497, 5.0))
        GeodesyXYZExt.fixbearing!(deg2rad(-90))

        enu = ENU(1.0, 0.0, 1.0)
        xyz = XYZ(enu)
        enu_ = ENU(xyz)
        @assert enu_ ≈ enu  # works :)

        lla = LLA(xyz)
        xyz_ = XYZ(lla)
        @assert xyz_ ≈ xyz  # works :)
    end
    @testset "block 3" begin
        @test_throws CoordinateSystemError XYZ(1., 2., 3.) + ENU(1., 2., 3.)
    end
    @testset "block 4" begin
        XYZ(1., 2., 3.) == ENU(1., 2., 3.)  # `true` (!)
    end
    @testset "block 5" begin
        using Unitful; using Unitful.DefaultSymbols
        Meters = typeof(1.0m)

        XYZ(1.0m, 2.0m, 3.0m) == XYZ(1., 2., 3.) * 1m
        ustrip.(uconvert.(mm, XYZ(1.0m, 2.0m, 3.0m))) == XYZ(1000., 2000., 3000.)
        @assert typeof(XYZ(1.0m, 2.0m, 3.0m)) == XYZ{Meters}
    end
end
