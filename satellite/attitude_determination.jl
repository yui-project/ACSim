using SatelliteToolbox
"""
sat_pos, sat_att = attitude_determination()

太陽センサと磁気センサから現在の位置・姿勢を決定

# Argments
 - `sun_vecs`：センサから受け取った太陽方向ベクトル
 - `mag_vecs`：センサから受け取った地磁場方向ベクトル
 - `gyr_vecs`：センサから受け取った角速度ベクトル

# Return
 - `sat_att`：衛星の推定姿勢クォータニオン
"""
function attitude_determination(sun_vecs::Vector,mag_vecs::Vector,gyr_vecs::Vector)
    sat_att=Quaternion(1,0,0,0)
    return sat_att
end
