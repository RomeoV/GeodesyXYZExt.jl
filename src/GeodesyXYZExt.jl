module GeodesyXYZExt
using Geodesy
import Geodesy: ENU
import Geodesy.CoordinateTransformations: Transformation
import Tau: τ
import Rotations: RotZ
import StaticArrays: FieldVector
import StaticArrays: map, _map, similar_type, Size

export ENUfromXYZ, XYZfromENU, XYZfromECEF, ECEFfromXYZ, XYZfromLLA, LLAfromXYZ,
       XYZfromUTM, UTMfromXYZ,
       XYZfromUTMZ, UTMZfromXYZ,
       XYZ

DATUM::Ref{Union{Datum, Nothing}} = nothing
ORIGIN::Ref{Union{LLA, Nothing}} = nothing
BEARING::Ref{Union{Float64, Nothing}} = nothing
function fixdatum!(datum::Union{Datum, Nothing})
    DATUM[] = datum
end
function fixorigin!(origin::Union{LLA, Nothing})
    ORIGIN[] = origin
end
function fixbearing!(bearing::Union{Float64, Nothing})
    BEARING[] = bearing
end

"""
    XYZ(x, y, z)

Alongtrack-Crosstrack-Elevation (XYZ) coordinates for the aviationo runway usecase. A local Cartesian coordinate system, linearized about a reference point.
"""
struct XYZ{T <: Number} <: FieldVector{3, T}
    x::T
    y::T
    z::T
end
# XYZ(x :: T, y :: T) where {T} = XYZ(x, y, zero(T))
@inline function XYZ(x,y,z)
    T = promote_type(promote_type(typeof(x),typeof(y)), typeof(z))
    XYZ{T}(x,y,z)
end
Base.show(io::IO, ::MIME"text/plain", xyz::XYZ) = print(io, "XYZ($(xyz.x), $(xyz.y), $(xyz.z))")



XYZ(xyz::XYZ, datum) = xyz
################################
### ENU <-> XYZ coordinates  ###
################################

ENU(xyz::XYZ) = ENU(xyz, BEARING[])
ENU(xyz::XYZ, bearing) = ENUfromXYZ(bearing)(xyz)
XYZ(enu::ENU) = XYZ(enu, BEARING[])
XYZ(enu::ENU, bearing) = XYZfromENU(bearing)(enu)

################################
### ECEF <-> XYZ coordinates ###
################################

XYZ(ecef::ECEF) = XYZ(ecef, ORIGIN[], BEARING[], DATUM[])
XYZ(ecef::ECEF, origin, bearing, datum) = XYZfromECEF(origin, bearing, datum)(ecef)
ECEF(xyz::XYZ) = ECEF(xyz, ORIGIN[], BEARING[], DATUM[])
ECEF(xyz::XYZ, origin, bearing, datum) = ECEFfromXYZ(origin, bearing, datum)(xyz)
################################
### LLA <-> XYZ coordinates ###
################################

XYZ(lla::LLA) =  XYZ(lla, ORIGIN[], BEARING[], DATUM[])
XYZ(lla::LLA, origin, bearing, datum) =  XYZfromLLA(origin, bearing, datum)(lla)
LLA(xyz::XYZ) = LLA(xyz, ORIGIN[], BEARING[], DATUM[])
LLA(xyz::XYZ, origin, bearing, datum) = LLAfromXYZ(origin, bearing, datum)(xyz)


################################
### XYZ <-> UTMZ coordinates ###
################################

XYZ(utm::UTMZ) = XYZ(utm, ORIGIN[], BEARING[], DATUM[])
XYZ(utm::UTMZ, origin, bearing, datum) = XYZfromUTMZ(origin, bearing, datum)(utm)
UTMZ(xyz::XYZ) = UTMZ(xyz, ORIGIN[], BEARING[], DATUM[])
UTMZ(xyz::XYZ, origin, bearing, datum) = UTMZfromXYZ(origin, bearing, datum)(xyz)



###############################
### XYZ <-> UTM coordinates ###
###############################

XYZ(utm::UTM, zone::Integer, hemisphere::Bool, origin, bearing, datum) = XYZfromUTM(origin, bearing, zone, hemisphere, datum)(utm)
UTM(xyz::XYZ, zone::Integer, hemisphere::Bool, origin, bearing, datum) = UTMfromXYZ(origin, bearing, zone, hemisphere, datum)(xyz)


##################
## ENU <-> XYZ ##
##################

"""
    ENUfromXYZ(xyz::XYZ)

Construct a `Transformation` object to convert from local `XYZ` coordinates
to local `ENU` coordinates centered at the same origin. This is a simple
permutation of coordinates and sign change for the altitude.

`bearing` is provided in rad, measured clockwise from north.
This follows the convention used in the runway description database.
"""
struct ENUfromXYZ <: Transformation
    bearing::Float64
end

Base.show(io::IO, ::ENUfromXYZ) = print(io, "ENUfromXYZ()")

function (obj::ENUfromXYZ)(xyz::XYZ)
    θ = τ/4 - obj.bearing  # counter-clockwise from E-axis
    R = RotZ(θ)
    convert(ENU, R*xyz)
end

"""
    XYZfromENU(bearing::Float64)

Construct a `Transformation` object to convert from local `ENU` coordinates
to local `XYZ` coordinates centered at the same origin. This is a simple
permutation of coordinates and sign change for the altitude.

See `ENUfromXYZ` for bearing documentation.
"""
struct XYZfromENU <: Transformation
    bearing::Float64
end     # singleton type

Base.show(io::IO, ::XYZfromENU) = print(io, "XYZfromENU()")

function (obj::XYZfromENU)(enu::ENU)
    θ = τ/4 - obj.bearing  # counter-clockwise from E-axis
    R = RotZ(θ)
    convert(XYZ, inv(R)*enu)
end

Base.inv(obj::ENUfromXYZ) = XYZfromENU(obj.bearing)
Base.inv(obj::XYZfromENU) = ENUfromXYZ(obj.bearing)



##################
## ECEF <-> XYZ ##
##################

"""
    XYZfromECEF(origin, bearing, datum)
    XYZfromECEF(origin::UTM, zone, isnorth, datum)
    XYZfromECEF(origin::ECEF, lat, lon)

Construct a composite transformation XYZfromENU(bearing) ∘ ENUfromECEF(origin, datum)
to convert from global `ECEF` coordinates to local `XYZ` coordinates centered at the `origin`.
This object pre-caches both the ECEF coordinates and latitude and longitude of the origin for maximal efficiency.
"""
XYZfromECEF(origin, bearing, datum) = XYZfromENU(bearing) ∘ ENUfromECEF(origin, datum)

"""
    ECEFfromXYZ(origin, bearing, datum)
    ECEFfromXYZ(origin::UTM, zone, isnorth, datum)
    ECEFfromXYZ(origin::ECEF, lat, lon)

Construct a composite transformation ECEFfromENU(origin,datum) ∘ ENUfromXYZ()
to convert from local `XYZ` coordinates centred at `origin` to global `ECEF` coodinates.
This object pre-caches both the ECEF coordinates and latitude and longitude of the origin for maximal efficiency.
"""
ECEFfromXYZ(origin, bearing, datum) = ECEFfromENU(origin, datum) ∘ ENUfromXYZ(bearing)



#################
## LLA <-> XYZ ##
#################

"""
    XYZfromLLA(origin, bearing, datum)

Creates composite transformation `XYZfromECEF(origin, bearing, datum) ∘ ECEFfromLLA(datum)`.
"""
XYZfromLLA(origin, bearing, datum) = XYZfromECEF(origin, bearing, datum) ∘ ECEFfromLLA(datum)

"""
    LLAfromXYZ(origin, bearing, datum)

Creates composite transformation `LLAfromECEF(datum) ∘ ECEFfromXYZ(origin, bearing, datum)`.
"""
LLAfromXYZ(origin, bearing, datum) = LLAfromECEF(datum) ∘ ECEFfromXYZ(origin, bearing, datum)




##################
## XYZ <-> UTMZ ##
##################

XYZfromECEF(origin::UTMZ, bearing, datum) = XYZfromECEF(LLAfromUTMZ(datum)(origin), datum)
ECEFfromXYZ(origin::UTMZ, bearing, datum) = ECEFfromXYZ(LLAfromUTMZ(datum)(origin), datum)

"""
    XYZfromUTMZ(origin, bearing, datum)

Creates composite transformation `ENUfromLLA(origin, bearing, datum) ∘ LLAfromUTMZ(datum)`.
"""
XYZfromUTMZ(origin, bearing, datum) = XYZfromLLA(origin, bearing, datum) ∘ LLAfromUTMZ(datum)

"""
    UTMZfromXYZ(origin, bearing, datum)

Creates composite transformation `UTMZfromLLA(datum) ∘ LLAfromXYZ(origin, bearing, datum)`.
"""
UTMZfromXYZ(origin, bearing, datum) = UTMZfromLLA(datum) ∘ LLAfromXYZ(origin, bearing, datum)

#################
## XYZ <-> UTM ##
#################
XYZfromECEF(origin::UTM, bearing, zone::Integer, isnorth::Bool, datum) = XYZfromECEF(LLAfromUTM(zone, isnorth, datum)(origin), bearing, datum)
ECEFfromXYZ(origin::UTM, bearing, zone::Integer, isnorth::Bool, datum) = ECEFfromXYZ(LLAfromUTM(zone, isnorth, datum)(origin), bearing, datum)

# Assume origin and utm point share the same zone and hemisphere
UTMfromXYZ(origin::UTM, bearing, zone::Integer, isnorth::Bool, datum) = UTMfromLLA(zone, isnorth, datum) ∘ LLAfromXYZ(UTMZ(origin, zone, isnorth), bearing, datum)
XYZfromUTM(origin::UTM, bearing, zone::Integer, isnorth::Bool, datum) = XYZfromLLA(UTMZ(origin, zone, isnorth), bearing, datum) ∘ LLAfromUTM(zone, isnorth, datum)

"""
    UTMfromXYZ(origin, bearing, zone, isnorth, datum)

Creates composite transformation `UTMfromLLA(zone, isnorth, datum) ∘ LLAfromXYZ(origin, bearing, datum)`.
If `origin` is a `UTM` point, then it is assumed it is in the given specified zone and hemisphere.
"""
UTMfromXYZ(origin, bearing, zone::Integer, isnorth::Bool, datum) = UTMfromLLA(zone, isnorth, datum) ∘ LLAfromXYZ(origin, bearing, datum)

"""
    XYZfromUTM(origin, bearing, zone, isnorth, datum)

Creates composite transformation `UTMfromLLA(zone, isnorth, datum) ∘ LLAfromXYZ(origin, bearing, datum)`.
If `origin` is a `UTM` point, then it is assumed it is in the given specified zone and hemisphere.
"""
XYZfromUTM(origin, bearing, zone::Integer, isnorth::Bool, datum) = XYZfromLLA(origin, bearing, datum) ∘ LLAfromUTM(zone, isnorth, datum)


### Base operations

# this basically also overloads (ENU(1., 2., 3.) * 1m)
# @inline function map(f, a1::ENU{<:Number})
#     ENU(_map(f, a1))
# end
# @inline function map(f, a1::XYZ{<:Number})
#     XYZ(_map(f, a1))
# end

similar_type(::Type{A}, ::Type{T}, s::Size{S}) where {A<:ENU, T, S} = ENU{T}
similar_type(::Type{A}, ::Type{T}, s::Size{S}) where {A<:XYZ, T, S} = XYZ{T}

end # module GeodesicXYZExt
