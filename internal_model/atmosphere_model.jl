using SatelliteToolbox
"""
x_ecef: Coordinate in the ECEF
"""
function atmospheric_models(x_ecef)
    x_geod = ECEFtoGeodetic(x_ecef)
    return expatmosphere(x_geod[3])
end
