include("sensors.jl")

"""
mag_vecs = [29.8522 -4.0303 35.8413]

mag_vols = mag_sensor(mag_vecs, "vol")
println(mag_vols)

"""

sun_vecs = [0.90 1.1111 0.383]
sun_pos = sun_sensor(sun_vecs, "z-")
println(sun_pos)
