using SatelliteToolbox
using LinearAlgebra

v = norm(V)//velosityvector norm

#require height of spacecraft
h = 
 398600 - 6378*v*v
#-----------------
/    v*v


#require atmosphere density
expatmosphere(h)