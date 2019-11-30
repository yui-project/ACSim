using SatelliteToolbox

#=*******************************
#noem of velocity vector v[km/s]
v = 7.66#norm(V)

#require height of spacecraft[km]
h = (398600 - 6378*v^2)/(v^2)
#conversion [km]->[m]
h = h*1e3
println(h)

#require atmosphere density[kg/m^3]
println(expatmosphere(h))
********************************=#

#v:norm of velocity vector v 
function atmospheric_models(v)
    h = (398600 - 6378*v^2)/(v^2)
    h = h*1e3
    return expatmosphere(h)
end
