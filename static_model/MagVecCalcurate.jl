include("Coordinate_Transform.jl")
using SatelliteToolbox

"""
@fn     MagVecCal
@input  date:計算したい時刻
        r_ecef : 地心座標系(ECEF)による位置3行1列[m]
@output 
"""
function MagVecCal(date::Number, r_ecef::AbstractVector)
    
    r_geodetic = ECEFtoGeodetic(r_ecef)
    # r_geocentric = GeodetictoGeocentric(r_geocentric)
    return igrf12(date, norm(r_ecef), r_geodetic[1], r_geodetic[2], show_warns = true)
end


r_ecef = [3.957729931663941e6,3.3091914769971482e6,3.737926200042794e6]
println(MagVecCal(2019, r_ecef))