include("sensors.jl")
using SatelliteToolbox

"""
q1 = SatelliteToolbox.Quaternion(cos(π/4), 1.0*sin(π/4), 0.0*sin(π/4), 0.0*sin(π/4))
q2 = SatelliteToolbox.Quaternion(cos(π/4), 0.0*sin(π/4), 0.0*sin(π/4), 1.0*sin(π/4))
dt = 5

ω = gyro_calcurate(q1, q2, dt)
println(ω)
"""

"""
mag_vecs = [29.8522 -4.0303 35.8413]

mag_vols = mag_sensor(mag_vecs, "vol")
println(mag_vols)

"""

sun_vecs = [-0.5 -0.5 1]
ss_out = sun_sensor(sun_vecs, "z-")
println(ss_out)


