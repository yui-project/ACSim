using Dates
using LinearAlgebra
using Quaternions

"""
sat_pos, sat_att = attitude_determination()

太陽センサと磁気センサから現在の位置・姿勢を決定

# Argments
 - `sun_vecs`：太陽方向ベクトル
 - `mag_vecs`：地磁場方向ベクトル

# Return
 - `sat_pos`：衛星の位置
 - `sat_att`：衛星の姿勢クォータニオン
"""
function attitude_determination()
    
end

"""
out_cur = attitude_control()

目標姿勢までPID制御する

# Argments
 - `sat_att`：現在の衛星姿勢クォータニオン
 - `sat_tar`：目標姿勢クォータニオン
 - `ω`：位置ベクトルの角速度

# Return
- `out_cur`：加える電流値
"""
function attitude_control()
    
end

"""
sat_pos, sat_att, out_cur = satellite()

姿勢決定、姿勢制御を行う

# Argments
 - `sun_vecs`：太陽方向ベクトル
 - `mag_vecs`：地磁場方向ベクトル
 - `sat_att`：現在の衛星姿勢クォータニオン
 - `sat_tar`：目標姿勢クォータニオン
 - `ω`：位置ベクトルの角速度

# Return
 - `sat_pos`：衛星の位置
 - `sat_att`：衛星の姿勢クォータニオン
"""
function satellite()
	
end