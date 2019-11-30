include("Coordinate_Transform.jl")
using SatelliteToolbox

"""
@fn     igrf
@input  date:計算したい時刻
        r : 地心座標系による位置3行1列[m]
@output 
"""
function igrf(date::Number, r::AbstractVector)
    
    r_geodetic = ECEFtoGeodetic(r)
    # r_geocentric = GeodetictoGeocentric(r_geocentric)
    return igrf12(date, norm(r), r_geodetic[1], r_geodetic[2], show_warns = true)
end


r = [3.957729931663941e6,3.3091914769971482e6,3.737926200042794e6]
println(igrf(2019, r))