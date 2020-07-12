include("orbit/orbit.jl")
include("external_model/external_model.jl")
include("external_model/disturbance_torque.jl")
include("internal_model/magnetic_torque.jl")
#include("dynamic_model/dynamic_model.jl")
include("dynamics/dynamics.jl")
#include("satellite/satellite.jl")
include("satellite/attitude_control.jl")
include("satellite/target_decision.jl")
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
	DataNum = 5000 #シミュレータ反復回数
	dt = 5 ##シミュレータの計算間隔 [s]
	start_time = DateTime(2019, 12, 19, 3, 27, 10)	#シミュレート開始時刻
	TLEFileName = "./orbit/ISS_TLE.txt"

	# 進行方向とSCOFx軸とのずれ（内積）
	dotvs = zeros(DataNum)

	# 軌道・太陽方向計算については先に行い、全 時間分を配列に保存する
	JD_log, x_ecef_log, v_ecef_log, x_geod_log, tles, eop_IAU2000A = orbit_cal(DataNum,dt,start_time,TLEFileName)
	sunvector_index = sunvector_model(JD_log[1], tles)
	#println(size(sunvector_index))

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
	disturbance = zeros(DataNum, 3)
	torqe = zeros(DataNum,3)
	mtq_currentlog = zeros(DataNum, 3)

	#衛星姿勢用変数
	sat_attqua_elements = zeros(DataNum+1, 4)
	sat_attqua_elements[1,:] = [cos(deg2rad(30)), 1/sqrt(3)*sin(deg2rad(30)), 1/sqrt(3)*sin(deg2rad(30)), 1/sqrt(3)*sin(deg2rad(30)),]
	sat_ω = zeros(DataNum+1, 3)
	sat_ω[1, :] = [0., 0., 0.]

	# トルク用変数
	airtorques = zeros(DataNum, 3)
	suntorques = zeros(DataNum, 3)
	magtorques = zeros(DataNum, 3)
	M_reqs = zeros(DataNum, 3)
	T_reqs = zeros(DataNum, 3)

	# 撮影用パラメータの設定
	limit_time = DataNum
	targetpos_geod = [rad2deg(-0.75), rad2deg(0.4), 25.7]
	targetpos_ecef = GeodetictoECEF(deg2rad(targetpos_geod[1]), deg2rad(targetpos_geod[2]), targetpos_geod[3])
	cam_viewangle = 40
	sat_axisval = 80
	cam_origindir = [0., 0., 1.]
	satqua = SatelliteToolbox.Quaternion(sat_attqua_elements[1,1], sat_attqua_elements[1,2], sat_attqua_elements[1,3], sat_attqua_elements[1,4])
	shoot_time, shoot_vec = shootingtime_decision2(x_ecef_log[1:DataNum, :], targetpos_ecef, 1, limit_time, cam_viewangle, sat_axisval)
	println("shoottime", shoot_time)
	println("shoot_vec", shoot_vec)
	targetqua, rotqua = targetqua_dicision(shoot_vec, satqua, sat_axisval)
	# targetqua = SatelliteToolbox.Quaternion(cos(deg2rad(5)), 0., sin(deg2rad(5)), 0.)
	println("targetqua:", targetqua)
	println("   rotqua:", rotqua)
	cam_dir = zeros(DataNum, 3)
	target_deffs = zeros(DataNum)

	for i=1:DataNum

		
		# 時刻表示
		current_time = start_time + Second(dt) * i
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
		mag_vec_seof = ecef_to_DCM(x_ecef_log[i,:],v_ecef_log[i,:],true)*mag_vec
		# mag_vec_seof = [0., 40000., 0.]
		mag_vecs[i,:] = mag_vec_seof
		println("         sun_vec:",sun_vec)
		sun_vecs[i,:] = sun_vec
		println("      atoms_dens:",atoms_dens)
		atoms_denses[i] = atoms_dens

		direct_on_SCOFs[i,:] = ecef_to_DCM(x_ecef_log[i,:],v_ecef_log[i,:],true) * v_ecef_log[i,:]


		# 衛星姿勢
		qua = SatelliteToolbox.Quaternion(sat_attqua_elements[i,1], sat_attqua_elements[i,2], sat_attqua_elements[i,3], sat_attqua_elements[i,4])
		

		
		#擾乱の計算
		v_scof = [norm(v_ecef_log[i,:]), 0., 0.]
		Ts = sun_pressure(sun_vecs[i,:], qua)
		Ts = [0., 0., 0.]
		suntorques[i,:] = Ts
		Ta = air_pressure(atoms_denses[i], v_scof, qua)
		Ta = [0., 0., 0.]
		airtorques[i,:] = Ta
		disturbance[i,:] = Ts + Ta

		#磁気トルカの計算
		
		#M = [0., 0., 0.]
		magvec_scsfqua = qua \ mag_vec_seof * qua
		atoms_denses[i] = atoms_dens
		# magvec_scsfqua = qua * mag_vecs[1,:] / qua
		magvec_scsf = [magvec_scsfqua.q1*10^(-9), magvec_scsfqua.q2*10^(-9), magvec_scsfqua.q3*10^(-9)]
		println("  scsf_mag_vec:", magvec_scsf)
		atoms_denses[i] = atoms_dens

		"""
		# B-dot法による回転抑制制御
		if i==1
			M = B_dot(magvec_scsf, sat_ω[i, :], sat_ω[i, :])
		else
			M = B_dot(magvec_scsf, sat_ω[i, :], sat_ω[i-1, :])
		end
		i_m = mm2current_theory(M)
		for j = 1:3
			if i_m[j] > 0.2
				i_m[j] = 0.2
			elseif i_m[j] < -0.2
				i_m[j] = -0.2
			end
		end
		println("  current_mag:",i_m)
		mtq_currentlog[i, :] = i_m 
		Tm = magnetic_torque(i_m, magvec_scsf)
		magtorques[i,:] = Tm
		"""

		# Cross-Product法による指向制御
		
		kp = 0.000000050
		kr = 0.0000050
		#=
		if i < DataNum/2
			tar_qua = targetqua2
		else
			tar_qua = targetqua
		end
		=#
		tar_qua = targetqua
		Treq, M = cross_product(tar_qua, qua, kp, kr, sat_ω[i, :], magvec_scsf)
		T_reqs[i, :] = Treq
		M_reqs[i, :] = M
		println("request_Moment:", M)
		
		i_m = mm2current_theory(M)
		for j = 1:3
			if i_m[j] > 0.2
				reduce_rate = 0.2 / i_m[j]
				i_m = reduce_rate * i_m
			elseif i_m[j] < -0.2
				reduce_rate = -0.2 / i_m[j]
				i_m = reduce_rate * i_m
			end
		end
		
		# M = [0.005, 0., 0.]
		# i_m = M
		println("  current_mag:",i_m)
		mtq_currentlog[i, :] = i_m
		Tm = magnetic_torque(i_m, magvec_scsf)
		magtorques[i,:] = Tm
		

		# ダイナミクス
		qua = SatelliteToolbox.Quaternion(sat_attqua_elements[i,1], sat_attqua_elements[i,2], sat_attqua_elements[i,3], sat_attqua_elements[i,4])
		I=[2*(0.1^2)/3 0. 0.;
		0. 2*(0.1^2)/3 0.;
		0. 0. 2*(0.1^2)/3]
		# next_qua, ω = dynamics(qua, sat_ω[i,:], Ta+Ts+Tm, I, dt)
		next_qua, ω = dynamics(qua, sat_ω[i,:], Treq, I, dt)
		sat_attqua_elements[i+1, :] = [next_qua.q0, next_qua.q1, next_qua.q2, next_qua.q3]
		sat_ω[i+1, :] = ω

		# カメラ方向
		cam_dir_q = qua * cam_origindir / qua
		cam_dir[i, :] = [cam_dir_q.q1, cam_dir_q.q2, cam_dir_q.q3]
		qua1 = qua * [0., 0., 1.] / qua
		vec1 = [qua1.q1, qua1.q2, qua1.q3]
		qua2 = targetqua * [0., 0., 1.] / targetqua
		vec2 = [qua2.q1, qua2.q2, qua2.q3]
		cθ = dot(vec1, vec2)
		sθ = norm(cross(vec1, vec2))
		target_deffs[i] = rad2deg(acos(cθ))
		if sθ < 0
			target_deffs[i] = target_deffs[i] * -1
		end
		

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
		plot_vec(disturbance, "disturbance")
		plot_3scalar(JD_log, sat_ω[1:DataNum, 1], sat_ω[1:DataNum, 2], sat_ω[1:DataNum, 3], ["x", "y", "z"], "omega")
		plot_3scalar(JD_log, airtorques[1:DataNum, 1], suntorques[1:DataNum, 1], magtorques[1:DataNum, 1], ["air", "sun", "mag"], "x_torques")
		plot_3scalar(JD_log, airtorques[1:DataNum, 2], suntorques[1:DataNum, 2], magtorques[1:DataNum, 2], ["air", "sun", "mag"], "y_torques")
		plot_3scalar(JD_log, airtorques[1:DataNum, 3], suntorques[1:DataNum, 3], magtorques[1:DataNum, 3], ["air", "sun", "mag"], "z_torques")
		plot_3scalar(JD_log, T_reqs[1:DataNum, 1], T_reqs[1:DataNum, 2], T_reqs[1:DataNum, 3], ["x", "y", "z"], "request_torques")
		plot_vec(T_reqs, "request_torques_vec")
		plot_3scalar(JD_log, mtq_currentlog[:, 1], mtq_currentlog[:, 2], mtq_currentlog[:, 3], ["x", "y", "z"], "MTQcurrent")
		plot_3scalar(JD_log, sat_attqua_elements[1:DataNum, 2], sat_attqua_elements[1:DataNum, 3], sat_attqua_elements[1:DataNum, 4], ["q1", "q2", "q3"], "attqua_1")
		plot_2scalar(JD_log, sat_attqua_elements[1:DataNum, 1], "attqua_2")
		plot_vec(cam_dir, "cam_dir")
		plot_2scalar(JD_log, target_deffs, "target_deffs")
		
		
	end

	
end

main()