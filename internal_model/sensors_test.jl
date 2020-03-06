include("sensors.jl")

"""
mag_vecs = [29.8522 -4.0303 35.8413]

mag_vols = mag_sensor(mag_vecs, "vol")
println(mag_vols)

"""
sun_vecs = [1.0 1.0 1.0]
sun_pos = sun_sensor(sun_vecs, "z-")
println(sun_pos)
