@testset "Co-ordinate system conversion" begin
    @testset "Fixed conversions" begin
        # ENU <-> XYZ
        @testset "bearing 90°" begin
            bearing = τ/4
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
        @testset "bearing 180°" begin
            bearing = τ/2
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

        @test XYZ(1., 2., 3.) + XYZ(1, 2, 3) == XYZ(2., 4., 6.)
        @test @fastmath XYZ(1., 2., 3.) + XYZ(1, 2, 3) == XYZ(2., 4., 6.)

        @test_throws CoordinateSystemError XYZ(1., 2., 3.) + ENU(1., 2., 3.)
        @test_throws CoordinateSystemError ENU(1., 2., 3.) + XYZ(1., 2., 3.)
        @test_throws CoordinateSystemError XYZ(1., 2., 3.) + ECEF(1., 2., 3.)
        @test_throws CoordinateSystemError @fastmath XYZ(1., 2., 3.) + ENU(1., 2., 3.)
        @test_throws CoordinateSystemError @fastmath ENU(1., 2., 3.) + XYZ(1., 2., 3.)
        @test_throws CoordinateSystemError @fastmath XYZ(1., 2., 3.) + ECEF(1., 2., 3.)
    end
end
