using SatelliteToolbox
"""
sat_pos, sat_att = attitude_determination()

太陽センサと磁気センサから現在の位置・姿勢を決定

# Argments
 - `(sun_vecs,mag_vecs,gyr_vecs)`：センサから受け取った太陽方向ベクトル,地磁場方向ベクトル,角速度ベクトルの組
 - `(sun_ina,mag_ina)`：現在位置での感性座標系での太陽方向ベクトルと地磁場方向ベクトルの組
 - `bef_quat`：1ステップ前の衛星の推定姿勢クォータニオン
# Return
 - `att_quat`：衛星の推定姿勢クォータニオン
"""
function attitude_determination((sun_vecs,mag_vec,gyr_vecs),(sun_ina,mag_ina),bef_quat)
    att_quat=Quaternion(1,0,0,0)
    return att_quat
end
