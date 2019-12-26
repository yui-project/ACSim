using SatelliteToolbox
using Plots
gr()

Lat = -π/2:10/180*π:π/2
x = zeros(18,3)
println(x)
anim = @animate for i=1:18
    println(i)
    x[i,:] = GeodetictoECEF(Lat[i],0,4000000)
    println(x[i,:])
    plot(x[1:i,1],x[1:i,3])
end
gif(anim, "test6_anim.gif", fps = 5)

#y = ECEFtoGeodetic([6367000, 0, -4000000])
#println(y)

#=
Lat = 30/180*π
Lon = 120/180*π
h = 4000000

sin_lat, cos_lat = sincos(Lat)
sin_lon, cos_lon = sincos(Lon)

a_wgs84  = 6378137.0
f_wgs84  = 1/298.257223563
b_wgs84  = a_wgs84*(1-f_wgs84)

N =a_wgs84/sqrt(1 - e_wgs84^2*sin_lat^2 )

x = [(N+h)*cos_lat*cos_lon, (N+h)*cos_lat*sin_lon, ((b_wgs84/a_wgs84)^2*N+h)*sin_lat]
println(x)
=#