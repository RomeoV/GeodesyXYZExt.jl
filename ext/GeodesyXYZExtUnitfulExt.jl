module GeodesyXYZExtUnitfulExt
import Geodesy: ENU
import GeodesyXYZExt: XYZ
import Unitful: ustrip, Units, Quantity

# ustrip(u::Units, a1::ENU{<:Quantity}) = map(x->ustrip(u, x), a1)
# ustrip(a1::ENU{<:Quantity}) = map(ustrip, a1)
# ustrip(u::Units, a1::XYZ{<:Quantity}) = map(x->ustrip(u, x), a1)
# ustrip(a1::XYZ{<:Quantity}) = map(ustrip, a1)
# ustrip(u::Units, a1::XYZ{<:Quantity}) = map(x->ustrip(u, x), a1)
# ustrip(a1::XYZ{<:Quantity}) = map(ustrip, a1)


end
