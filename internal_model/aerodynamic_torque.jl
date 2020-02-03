using SatelliteToolbox
using Plots
"""
available data are height, attitude, velocity, center of aeroforce
"""
function aerodynamics_torque(h, a, v, ca)