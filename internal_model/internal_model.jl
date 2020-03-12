include("aerodynamic_torque.jl")
include("magnetic_torque.jl")
include("sensors.jl")

using LinearAlgebra
using SatelliteToolbox
#=引数候補
i_m = [1, 1, 1]
B = [1, 1, 1]
rho = 1
v = 1
e_r = [1, 0, 0]
sun_vecs = [-0.5 -0.5 1.]
v_ecef_pu = [1. 1. 1.; 1. 1. 1.]
attqua1 = SatelliteToolbox.Quaternion(1.0/sqrt(2.0), 1.0/sqrt(2.0), 0.0, 0.0)
attqua2 = SatelliteToolbox.Quaternion(1.0/sqrt(2.0), 0.0, 1.0/sqrt(2.0), 0.0)
dt = 5
=#
"""
torque, ss_out, ms_out, acce_out, gyro_out = internal_model(i_m, B, rho, v, e_r, sun_vecs, v_ecef_pu, attqua1, attqua2, dt)

# Argments
- `i_m`：各磁気トルカに流れる電流
- `B`：磁束密度ベクトル＠SCSF
- `rho`：大気密度
- `v`：大気に対する衛星の速度
- `e_r`：衛星の進行方向の単位ベクトル@SCSF（body）
- `sun_vecs`：太陽方向ベクトル@SCSF (body)
- `v_ecef_pu`：今回ループと1回前ループにおける衛星速度ベクトル@ECEF (3*2 1行目が1回前、2行目が今回)
- `attqua1`：2回前ループにおける衛星姿勢クォータニオン
- `attqua2`：1回前ループにおける衛星姿勢クォータニオン
- `dt`：ループの時間間隔

# Returns
- `torqe`：衛星にかかるトルク@SCSF
- `ss_out`：太陽センサが出力する電流値
- `ms_out`：磁気センサが出力する磁束密度データ（I2Cシリアル出力想定）
- `acce_out`：ジャイロセンサが出力する加速度データ@ECEF（I2Cシリアル出力想定）
- `gyro_out`：ジャイロセンサが出力する角速度データ@SCSF（I2Cシリアル出力想定）
"""
function internal_model(i_m, B, rho, v, e_r, sun_vecs, v_ecef_pu, attqua1, attqua2, dt)
	#=
	設定パラメータ
	=#
	C_d = 1 #衛星のサーフェースの抵抗係数
	sur_num = 6 #衛星の面数
	p_c = [150 0 0; -150 0 0; 0 150 0; 0 -150 0; 0 0 150; 0 0 -150] #各面の圧力中心
	e_a = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1] #各面の法線ベクトル
	a = [90000, 90000, 90000, 90000, 90000, 90000] #衛星の各面積
	torque = zeros(3)

	#トルク計算
	aer_tor = aerodynamic_torque(sur_num, C_d, rho, v, a, e_r, e_a, p_c)
	mag_tor = magnetic_torque(i_m, B)
	torque = aer_tor + mag_tor

	#センサ出力
	ss_out = sun_sensor(sun_vecs, "z-")
	ms_out = mag_sensor(B)
	acce_out, gyro_out = gyro_sensor(v_ecef_pu, attqua1, attqua2, dt)
	
	return torque, ss_out, ms_out, acce_out, gyro_out
end