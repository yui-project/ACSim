include("aerodynamic_torque.jl")
include("magnetic_torque.jl")

using SatelliteToolbox

#=
設定パラメータ
=#
C_d = 1 #衛星のサーフェースの抵抗係数
sur_num = 6 #衛星の面数
p_c = [150 0 0; -150 0 0; 0 150 0; 0 -150 0; 0 0 150; 0 0 -150] #各面の圧力中心
e_a = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1] #各面の法線ベクトル
a = [90000, 90000, 90000, 90000, 90000, 90000] #衛星の各面積

#引数候補
i_m = [1, 1, 1]
B = [1, 1, 1]
rho = 1
v = 1
e_r = [1, 0, 0] 
"""
torque, ??? = internal_model(current_time, x_ecef, v_ecef, x_geod, ...)

# Argments
- ''：
- ''：
- ''：
- ''：
- ''：
- ''：

# Returns
- 'torqe'：衛星にかかるトルク@SCSF
- '???'：
"""
function internal_model(current_time, x_ecef, v_ecef, x_geod)
	torque = zeros(3)
	aer_tor = aerodynamic_torque(C_d, rho, v, a, e_r, e_a, p_c)
	mag_tor = magnetic_torque(i_m, B)
	torque = aer_tor + mag_tor
	return torque
end