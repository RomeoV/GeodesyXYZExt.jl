import Unitful.DefaultSymbols: m, °, rad, mm
import Unitful: ustrip, uconvert
Angle = Union{typeof(1.0°), typeof(1.0rad)};
Meters = typeof(1.0m)

@testset "unitful" begin
    @testset "cast unitless to unitful" begin
        @test typeof(XYZ(1., 2., 3.) * 1m) == XYZ{Meters}
        @test typeof(ENU(1., 2., 3.) * 1m) == ENU{Meters}
    end
    @testset "ustrip" begin
        @test typeof(ustrip.(XYZ(1., 2., 3.) * 1m)) == XYZ{Float64}
        @test typeof(ustrip.(ENU(1., 2., 3.) * 1m)) == ENU{Float64}
        @test ustrip.(m, ENU(1., 2., 3.) * 1m) == ENU{Float64}(1., 2., 3.)
        @test ustrip.(m, XYZ(1., 2., 3.) * 1m) == XYZ{Float64}(1., 2., 3.)
        @test ustrip.(mm, ENU(1., 2., 3.) * 1m) == ENU{Float64}(1., 2., 3.)*1000
        @test ustrip.(mm, XYZ(1., 2., 3.) * 1m) == XYZ{Float64}(1., 2., 3.)*1000
        @test uconvert.(mm, XYZ(1., 2., 3.) * 1m) == XYZ(1.0mm, 2.0mm, 3.0mm)*1000
    end
end
