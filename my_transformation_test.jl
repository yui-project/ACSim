include("my_transformation.jl")
using SatelliteToolbox

"""
test for eceftoseof
"""
#=
satpos = [6771000 0 0]
satvel = [0 8000 0]

vin = [0. 0. 0.]
vout = eceftoseof(vin, satpos, satvel, "pos")
println(vout)
=#

"""
test for seoftoscsf
"""

vin = [1. 0. 0.]
qua = SatelliteToolbox.Quaternion(cos(π/4), sin(π/4), 0., 0.)
vout = seoftoscsf(vin, qua)
println(vout)
