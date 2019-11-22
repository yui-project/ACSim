using SatelliteToolbox

#velosityvector norm[km/s]
v = 7.66#norm(V)

#require height of spacecraft[km]
h = (398600 - 6378*v^2)/(v^2)
#conversion [km]->[m]
h = h*1e3
println(h)

#require atmosphere density[kg/m^3]
println(expatmosphere(h))

