include("orbit/orbit.jl")
include("external_model/external_model.jl")
#include("dynamic_model/dynamic_model.jl")
#include("dynamics/dynamics.jl")
#include("satellite/satellite.jl")

using SatelliteToolbox
using Dates


#=
設定パラメータ
=#
DataNum =1000 #シミュレータ反復回数
dt = 5 ##シミュレータの計算間隔 [s]
start_time = DateTime(2019, 12, 19, 3, 27, 10)	#シミュレート開始時刻
TLEFileName = "./orbit/ISS_TLE.txt"

# 軌道計算については先に行い、全時間分を配列に保存する
JD_log, x_ecef_log , x_geod_log, v_ecef_log= orbit_cal(DataNum,dt,start_time,TLEFileName)


for i=1:10#DataNum
	current_time = start_time+Second(dt) * i
	println("\nCurrent DateTime:",current_time)
	println("        JD:",JD_log[i])
	println("x_ecef_log:",x_ecef_log[i,:])

	# 衛星外環境モデル
	mag_vel, sun_vec, atoms_dens = external_model(current_time,x_ecef_log[i,:])
	println("mag_vel:",mag_vel)
	println("sun_vec:",sun_vec)
	println("atoms_dens：",atoms_dens)


end