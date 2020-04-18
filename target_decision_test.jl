include("orbit/orbit.jl")
#include("external_model/external_model.jl")
#include("dynamic_model/dynamic_model.jl")
#include("dynamics/dynamics.jl")
#include("satellite/satellite.jl")
#include("plot/plot_plots.jl")
#include("plot/plot_makie.jl")
include("coordinates.jl")
include("satellite/target_decision.jl")

using SatelliteToolbox
using Dates
using LinearAlgebra
using Plots
gr()

using HDF5

function main()
	#=
	設定パラメータ
	=#
	DataNum =1500 #シミュレータ反復回数
	dt = 30 ##シミュレータの計算間隔 [s]
	start_time = DateTime(2020, 04, 09, 15, 18, 00)	#シミュレート開始時刻
	TLEFileName = "./orbit/ISS_TLE.txt"

	# 進行方向とSCOFx軸とのずれ（内積）
	dotvs = zeros(DataNum)

	# 撮影地点
	target_geod = [36.0834 140.0766 0.]
	limit_time = 1500
	print("target_geod : ")
	print(deg2rad(target_geod[1]))
	print(", ")
	print(deg2rad(target_geod[2]))
	print(", ")
	println(deg2rad(target_geod[3]))

	# カメラパラメータ
	cam_viewangle = 40
	sat_axisval = 40

	#=
	# 軌道計算については先に行い、全  時間分を配列に保存する
	JD_log, x_ecef_log, v_ecef_log, x_geod_log = orbit_cal(DataNum,dt,start_time,TLEFileName)

	h5open("orbit.h5", "w") do file
		write(file, "JD_log", JD_log)  # alternatively, say "@write file A"
		write(file, "x_ecef_log", x_ecef_log)
		write(file, "v_ecef_log", v_ecef_log)
		write(file, "x_geod_log", x_geod_log)
	end
	=#
	JD_log = h5read("orbit.h5", "JD_log")
	x_ecef_log = h5read("orbit.h5", "x_ecef_log")
	v_ecef_log = h5read("orbit.h5", "v_ecef_log")
	x_geod_log = h5read("orbit.h5", "x_geod_log")
	

	println(size(x_ecef_log))
	target_geod = [deg2rad(target_geod[1]), deg2rad(target_geod[2]), target_geod[3]]
	target_ecef = GeodetictoECEF(target_geod[1], target_geod[2], target_geod[3])

	currenttime = 1
	shootingtime, shooting_vec = shootingtime_decision2(x_ecef_log[currenttime:currenttime+limit_time-1, :], target_ecef, currenttime, limit_time, cam_viewangle, sat_axisval)
	println(shootingtime)
	println(shooting_vec)

	M = ecef_to_DCM(x_ecef_log[shootingtime,:],v_ecef_log[shootingtime,:])
	shooting_vec_seof = M*shooting_vec
	
	current_attqua = SatelliteToolbox.Quaternion(1., 0., 0., 0.)
	targetqua, rotqua = targetqua_dicision(shooting_vec_seof, current_attqua, sat_axisval)
	println(targetqua)
	println(rotqua)
	
	# 描画のための処理
	target_geod_plot = [target_geod[1] target_geod[2]; target_geod[1] target_geod[2]]
	shoot_geod_plot = [x_geod_log[shootingtime,1] x_geod_log[shootingtime,2]; x_geod_log[shootingtime,1] x_geod_log[shootingtime,2]]

	
	cθ = dot(shooting_vec, -x_ecef_log[shootingtime, :])/(norm(shooting_vec)*norm(x_ecef_log[shootingtime, :]))
	tθ = sqrt(1/(cθ^2)-1)
	ta = tand(cam_viewangle/2)
	tθa = (tθ+ta)/(1-tθ*ta)
	shootcentor_xdis = shooting_vec[1]*shooting_vec[3]/x_geod_log[shootingtime, 3]
	shootcentor_ydis = shooting_vec[2]*shooting_vec[3]/x_geod_log[shootingtime, 3]
	shootcentor = ECEFtoGeodetic(x_ecef_log[shootingtime, :] + [shootcentor_xdis, shootcentor_ydis, 0.])
	shootrange = tθa*x_geod_log[shootingtime, 3]
	shootendpoint = ECEFtoGeodetic(x_ecef_log[shootingtime, :] - [shootrange, 0., x_ecef_log[shootingtime, 3]])
	Δx = shootendpoint[1] - shootcentor[1]
	Δy = shootendpoint[2] - shootcentor[2]
	r = sqrt(Δx^2 + Δy^2)
	shootedge = zeros(37,2)
	for i = 1:37
		shootedge[i, :] = [shootcentor[2] + r*cosd(i*10), shootcentor[1] + r*sind(i*10)]
	end
	

	# 描画
	plot(x_geod_log[:,2], x_geod_log[:, 1], legend=:topleft)
	plot!(target_geod_plot[:,2], target_geod_plot[:,1], seriestype=:scatter, legend=:topleft)
	plot!(shoot_geod_plot[:, 2], shoot_geod_plot[:,1], seriestype=:scatter, legend=:topleft)
	plot!(shootedge[:, 1], shootedge[:, 2], legend=:topleft)	
	#plot!(target_ecef_plot[:, 1], target_ecef_plot[:, 2], target_ecef_plot[:, 3], seriestype=:scatter)
	xlabel!("x")
	ylabel!("y")
	#zlabel!("z")

	savefig("./figs/TargetDecisionTest2.png")

end

main()