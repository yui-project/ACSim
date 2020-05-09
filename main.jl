include("orbit/orbit.jl")
include("external_model/external_model.jl")
#include("dynamic_model/dynamic_model.jl")
#include("dynamics/dynamics.jl")
#include("satellite/satellite.jl")
include("plot/plot_plots.jl")
#include("plot/plot_makie.jl")
include("coordinates.jl")

using SatelliteToolbox
using Dates
using LinearAlgebra

using HDF5

function main()
	#=
	設定パラメータ
	=#
	DataNum =1000 #シミュレータ反復回数
	dt = 5 ##シミュレータの計算間隔 [s]
	start_time = DateTime(2019, 12, 19, 3, 27, 10)	#シミュレート開始時刻
	TLEFileName = "./orbit/ISS_TLE.txt"

	# 進行方向とSCOFx軸とのずれ（内積）
	dotvs = zeros(DataNum)

	# 軌道・太陽方向計算については先に行い、全 時間分を配列に保存する
	JD_log, x_ecef_log, v_ecef_log, x_geod_log, tles, eop_IAU2000A = orbit_cal(DataNum,dt,start_time,TLEFileName)
	sunvector_index = sunvector_model(JD_log[1], tles)

	h5open("orbit.h5", "w") do file
		write(file, "JD_log", JD_log)  # alternatively, say "@write file A"
		write(file, "x_ecef_log", x_ecef_log)
		write(file, "v_ecef_log", v_ecef_log)
		write(file, "x_geod_log", x_geod_log)
		write(file, "sunvector_index", sunvector_index)
	end

	# 初期状態の衛星位置@ECI
	x_eci_1st = rECEFtoECI(ITRF(), GCRF(), JD_log[1], eop_IAU2000A)*x_ecef_log[1,:]

	# SCOF上での衛星の進行方向
	direct_on_SCOFs = zeros(DataNum,3)

	# 衛星外環境モデル用変数
	mag_vecs = zeros(DataNum,3)
	sun_vecs = zeros(DataNum,3)
	atoms_denses = zeros(DataNum)

  	# 衛星内環境モデル用変数
  	torqe = zeros(DataNum,3)

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
		mag_vec, sun_vec, atoms_dens = external_model(current_time,x_ecef_log[i,:],x_geod_log[i,:], x_eci_1st, eop_IAU2000A, sunvector_index)
		println("         mag_vec:",mag_vec)
		mag_vecs[i,:] = ecef_to_DCM(x_ecef_log[i,:],v_ecef_log[i,:],true)*mag_vec
		println("         sun_vec:",sun_vec)
		sun_vecs[i,:] = sun_vec
		println("      atoms_dens:",atoms_dens)
		atoms_denses[i] = atoms_dens

		direct_on_SCOFs[i,:] = ecef_to_DCM(x_ecef_log[i,:],v_ecef_log[i,:],true) * v_ecef_log[i,:]



	end

	plot = true

	if plot == true
		plot_2scalar(JD_log,dotvs,"dotvs")
		plot_vec(direct_on_SCOFs,"direct_on_SCOFs")
		
		plot_2scalar(JD_log,atoms_denses,"atoms_dens")
		plot_2scalar(JD_log,x_geod_log[:,3],"height")
		plot_vec(mag_vecs,"mag_vec")
		plot_vec(sun_vecs,"sun_vec")
		plot_vec(x_ecef_log,"x_ecef_log")
		plot_vec(v_ecef_log,"v_ecef_log")
		plot_vec(x_geod_log,"x_geod_log")
		plot_2scalar(x_geod_log[:,2],x_geod_log[:,1],"x_geod_log_2d")
	end

	


	
end

main()