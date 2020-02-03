include("orbit/orbit.jl")
include("external_model/external_model.jl")
#include("dynamic_model/dynamic_model.jl")
#include("dynamics/dynamics.jl")
#include("satellite/satellite.jl")
include("graph_plot.jl")

using SatelliteToolbox
using Dates
using LinearAlgebra

function main()
	#=
	設定パラメータ
	=#
	DataNum =1000 #シミュレータ反復回数
	dt = 5 ##シミュレータの計算間隔 [s]
	start_time = DateTime(2019, 12, 19, 3, 27, 10)	#シミュレート開始時刻
	TLEFileName = "./orbit/ISS_TLE.txt"

	# 軌道計算については先に行い、全時間分を配列に保存する
	JD_log, x_ecef_log, v_ecef_log, x_geod_log = orbit_cal(DataNum,dt,start_time,TLEFileName)
	dotvs = zeros(DataNum)

	# 衛星外環境モデル用変数
	mag_vecs = zeros(DataNum,3)
	sun_vecs = zeros(DataNum,3)
	atoms_denses = zeros(DataNum)

	for i=1:DataNum

		# 時刻表示
		current_time = start_time+Second(dt) * i
		println("")
		println("Current DateTime:",current_time)
		println("              JD:",JD_log[i])
		println("      x_ecef_log:",x_ecef_log[i,:])
		println("      v_ecef_log:",v_ecef_log[i,:])
		println("      x_geod_log:",x_geod_log[i,:])
		dotv = dot(x_ecef_log[i,:] , v_ecef_log[i,:])
		dotvs[i] = dotv
		println("            dotv:",dotv)

		# 衛星外環境モデル計算
		mag_vec, sun_vec, atoms_dens = external_model(current_time,x_ecef_log[i,:],x_geod_log[i,:])
		println("         mag_vec:",mag_vec)
		mag_vecs[i,:] = mag_vec
		println("         sun_vec:",sun_vec)
		sun_vecs[i,:] = sun_vec
		println("      atoms_dens:",atoms_dens)
		atoms_denses[i] = atoms_dens


		

	end

	plot_2scalar(JD_log,dotvs,"dotvs")
	#=
	plot_2scalar(JD_log,atoms_denses,"atoms_dens")
	plot_2scalar(JD_log,x_geod_log[:,3],"height")
	plot_vec(mag_vecs,"mag_vec")
	plot_vec(sun_vecs,"sun_vec")
	plot_vec(x_ecef_log,"x_ecef_log")
	plot_vec(v_ecef_log,"v_ecef_log")
	plot_vec(x_geod_log,"x_geod_log")
	plot_2scalar(x_geod_log[:,2],x_geod_log[:,1],"x_geod_log_2d")
	=#


	
end

main()