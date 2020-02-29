include("aerodynamic_torque.jl")
include("magnetic_torque.jl")
include("sensor.jl")

using LinearAlgebra
#=引数候補
i_m = [1, 1, 1]
B = [1, 1, 1]
rho = 1
v = 1
e_r = [1, 0, 0]
=#
"""
torque = internal_model(i_m, B, rho, v, e_r)

# Argments
- `i_m`：各磁気トルカに流れる電流
- `B`：磁束密度ベクトル＠SCSF
- `rho`：大気密度
- `v`：大気に対する衛星の速度
- `e_r`：衛星の進行方向の単位ベクトル@SCSF（body）
- ``：

# Returns
- `torqe`：衛星にかかるトルク@SCSF
- `sun_vol`：太陽センサが出力する電圧
- `mag_vol`：磁気センサが出力する電圧
"""
function internal_model(i_m, B, rho, v, e_r)
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
	sun_vol = sun_sensor(sun_vecs, sat_pos, sat_att)
	mag_vol = mag_sensor(mag_vecs, sat_pos, sat_att)
	
	return torque, sun_vol, mag_vol
end