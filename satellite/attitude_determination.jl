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


"""
satqua, satvec, omega = attitute_empty(satqua, satvec, omega)

衛星姿勢をそのまま返す空関数（+ 誤差を与える）

# Argments
 - `satqua` : 衛星姿勢
 - `satvec` : 衛星速度@ECEF
 - `omega`  : 衛星回転速度

"""
function attitude_empty(satqua, satvec, omega)

    return satqua, satvec, omega
end

"""
roll, pitch, yaw = EulerAngles(quaternion)

クォータニオンを，X,Y,Z軸の順で回転させた際のオイラー角に変換する関数

# Argments
 - `quaternion` : クォータニオン
# Return
 - `roll`  : ロール角
 - `pitch` : ピッチ角
 - `yaw`   : ヨー角
 """
function  Euler_Angles(qua)
    q0 = qua.q0
    q1 = qua.q1
    q2 = qua.q2
    q3 = qua.q3

    roll = atan(2*(q0*q1+q2*q3)/(q0^2-q1^2-q2^2+q3^2))
    # println(typeof(roll))
    pitch = asin(2*(q0*q2-q1*q3))
    yaw = atan(2*(q0*q3+q1*q2)/(q0^2+q1^2-q2^2-q3^2))

    return roll, pitch, yaw
end 