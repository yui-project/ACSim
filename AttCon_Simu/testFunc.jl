a_wgs84=6378137.
f_wgs84 = 1/298.257223563
b_wgs84 = a_wgs84*(1-f_wgs84)
e_wgs84 = sqrt((a_wgs84^2-b_wgs84^2)/a_wgs84^2)

#export ECEFtoGeodetic, GeodetictoECEF
export GeodetictoECEF

function GeodetictoECEF2(lat::Number, lon::Number, h::Number)
    # Auxiliary variables.
    sin_lat, cos_lat = sincos(lat)
    sin_lon, cos_lon = sincos(lon)

    # Radius of curvature [m].
    N = a_wgs84/sqrt(1 - e_wgs84^2*sin_lat^2 )

    # Compute the position in ECEF frame.
    SVector{3}((                     N + h)*cos_lat*cos_lon,
               (                     N + h)*cos_lat*sin_lon,
               ( (b_wgs84/a_wgs84)^2*N + h)*sin_lat)
    
    println(SVector)
end