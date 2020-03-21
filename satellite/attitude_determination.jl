using SatelliteToolbox
using LinearAlgebra
"""
sat_pos, sat_att = attitude_determination()

太陽センサと磁気センサから現在の位置・姿勢を決定

# Argments
 - `(sun_vecs,mag_vecs,gyr_vecs)`：センサから受け取った太陽方向ベクトル,地磁場方向ベクトル,角速度ベクトルの組
 - `(isun_vecs,imag_vecs)`：現在位置での慣性座標系での太陽方向ベクトルと地磁場方向ベクトルの組
 - `bef_quat`：1ステップ前の衛星の推定姿勢クォータニオン
 -`dt`：刻み幅
# Return
 - `att_quat`：衛星の推定姿勢クォータニオン
"""
function attitude_determination((sun_vecs,mag_vec,gyr_vecs),(isun_vecs,imag_vecs),bef_quat,dt)
    att_quat=bef_quat+0.5bef_quat*gyr_vecs*dt
    att_quat/=norm(att_quat)
    return att_quat
end
