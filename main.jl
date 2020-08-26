include("orbit/orbit.jl")
include("external_model/external_model.jl")
include("external_model/disturbance_torque.jl")
include("internal_model/magnetic_torque.jl")
#include("dynamic_model/dynamic_model.jl")
include("dynamics/dynamics.jl")
#include("satellite/satellite.jl")
include("satellite/attitude_control.jl")
include("satellite/attitude_determination.jl")
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
	sat_attqua_elements[1,:] = [cos(deg2rad(35)), -sin(deg2rad(35))/sqrt(2), sin(deg2rad(35))/sqrt(2), 0.]
	sat_ω = zeros(DataNum+1, 3)
	sat_ω[1, :] = [0., 0., 0.]

	# トルク用変数
	airtorques = zeros(DataNum, 3)
	suntorques = zeros(DataNum, 3)
	magtorques = zeros(DataNum, 3)
	M_reqs = zeros(DataNum, 3)
	T_reqs = zeros(DataNum, 3)
	Tms = zeros(DataNum, 3)
	Terrors = zeros(DataNum, 3)

	# 誤差検証
	x_ecef_error = [0., 0., 0. ]
	targetqua_error = SatelliteToolbox.Quaternion(cos(deg2rad(5)), 0, sin(deg2rad(5)), 0)
	

	# 撮影用パラメータの設定
	limit_time = DataNum
	targetpos_geod = [0., -50., 25.7]
	targetpos_ecef = GeodetictoECEF(deg2rad(targetpos_geod[1]), deg2rad(targetpos_geod[2]), targetpos_geod[3])
	cam_viewangle = 40
	sat_axisval = 80
	cam_origindir = [0., 0., 1.]
	satqua = SatelliteToolbox.Quaternion(sat_attqua_elements[1,1], sat_attqua_elements[1,2], sat_attqua_elements[1,3], sat_attqua_elements[1,4])
	shoot_time = shootingtime_decision2(x_ecef_log[1:DataNum, :], targetpos_ecef, 1, limit_time, cam_viewangle, sat_axisval)
	println("shoottime", shoot_time)
	targetqua = targetqua_dicision(targetpos_ecef, x_ecef_log[shoot_time-120, :], v_ecef_log[shoot_time-120, :], sat_axisval)
	# targetqua = SatelliteToolbox.Quaternion(cos(deg2rad(44)), sin(deg2rad(44))/sqrt(2), sin(deg2rad(44))/sqrt(2), 0.)
	println("targetqua:", targetqua)
	cam_dir = zeros(DataNum, 3)
	tar_dir = zeros(DataNum, 3)
	target_deffs = zeros(DataNum)
	targetlocus_picture = zeros(240, 2)
	campos_log = zeros(239, 2)
	tarpos_log = zeros(239, 2)
	j = 1

	controlmethod = zeros(DataNum)

	# sat_attqua_elements[1, : ] = [targetqua.q0, targetqua.q1, targetqua.q2, targetqua.q3]


	
	for i=1:DataNum

		# 時刻表示
		current_time = start_time + Second(dt) * i
		println("")
		println("Current DateTime:",current_time)
		println("       loop time:", i)
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
		# Ts = [0., 0., 0.]
		suntorques[i,:] = Ts
		Ta = air_pressure(atoms_denses[i], v_scof, qua)
		# Ta = [0., 0., 0.]
		airtorques[i,:] = Ta
		disturbance[i,:] = Ts + Ta

		#M = [0., 0., 0.]
		magvec_scsfqua = qua \ mag_vec_seof * qua
		atoms_denses[i] = atoms_dens
		# magvec_scsfqua = qua * mag_vecs[1,:] / qua
		magvec_scsf = [magvec_scsfqua.q1*10^(-9), magvec_scsfqua.q2*10^(-9), magvec_scsfqua.q3*10^(-9)]
		println("  scsf_mag_vec:", magvec_scsf)
		atoms_denses[i] = atoms_dens


		# 姿勢決定による衛星姿勢情報の取得
		qua_ad, vel_ad, omega_ad = attitude_empty(qua, v_ecef_log[i, :], sat_ω[i, :])

		#=
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
		=#

		# Cross-Product法による指向制御
		#=
		kp = 0.000030
		kr = 0.00030
		#=
		if i < DataNum/2
			tar_qua = targetqua2
		else
			tar_qua = targetqua
		end
		=#
		tar_qua = targetqua
		Treq, M = cross_product(tar_qua, qua_ad, kp, kr, omega_ad, magvec_scsf)
		# T_reqs[i, :] = Treq
		=#

		kp = 0.000030
		kr = 0.00030
		
		#=
		tar_qua = targetqua
		if i == 1 
			M = B_dot(magvec_scsf, sat_ω[i, :], sat_ω[i, :])
			controlmethod[i] = -1
			Treq = cross(M, magvec_scsf)
			Terror = [0., 0., 0.]

		elseif norm(shoot_time - i) > 150  #&& norm(shoot_time-540 - i) > 150
			M = B_dot(magvec_scsf, sat_ω[i, :], sat_ω[i-1, :])
			controlmethod[i] = -1
			Treq = cross(M, magvec_scsf)
			Terror = [0., 0., 0.]
		
		else
			M = B_dot(magvec_scsf, sat_ω[i, :], sat_ω[i-1, :])
			controlmethod[i] = 1
			Treq = cross(M, magvec_scsf)
			Terror = [0., 0., 0.]
		
		end
		Terrors[i, :] = Terror
		=#
		
		tarqua = targetqua_error * targetqua
		Treq, M = cross_product(tarqua,qua_ad, kp, kr, omega_ad, magvec_scsf)
		# Treq, M, n = crossproduct_adj(targetqua, qua_ad, kp, kr, omega_ad, magvec_scsf)
		# Terrors[i, :] = dot(Treq, magvec_scsf)/(norm(magvec_scsf)^2) * magvec_scsf
		# controlmethod[i] = n

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
		
		# Treqの離散化
		Tmax = 1.0*10^(-5) # 最大出力トルク
		t_scatternum = 255 # Treqを±255段階に分ける
		Tresol = Tmax / t_scatternum # 出力トルクの分解能

		Treq = round.(Treq / Tresol) * Tresol
		T_reqs[i, :] = Treq


		
		

		# ダイナミクス
		qua = SatelliteToolbox.Quaternion(sat_attqua_elements[i,1], sat_attqua_elements[i,2], sat_attqua_elements[i,3], sat_attqua_elements[i,4])
		I=[(0.1^2)/6 0.        0.;
		    0.       (0.1^2)/6 0.;
		    0.       0.        (0.1^2)/6]
		# next_qua, ω = dynamics(qua, sat_ω[i,:], Ta+Ts+Tm, I, dt)
		# next_qua, ω = dynamics(qua, sat_ω[i,:], Tm, I, dt)
		# next_qua, ω = dynamics(qua, sat_ω[i,:], Treq, I, dt)
		next_qua, ω = dynamics(qua, sat_ω[i,:], Treq+Ta+Ts, I, dt)
		sat_attqua_elements[i+1, :] = [next_qua.q0, next_qua.q1, next_qua.q2, next_qua.q3]
		sat_ω[i+1, :] = ω
		println("quaternion:", qua)
		println("onega : ", ω)

		# カメラ方向
		cam_dir_q = qua * cam_origindir / qua
		cam_dir[i, :] = [cam_dir_q.q1, cam_dir_q.q2, cam_dir_q.q3]
		qua1 = qua * cam_origindir / qua
		vec1 = [qua1.q1, qua1.q2, qua1.q3]
		qua2 = targetqua * cam_origindir / targetqua
		vec2 = [qua2.q1, qua2.q2, qua2.q3]
		cθ = dot(vec1, vec2) / norm(vec1) / norm(vec2)
		if cθ > 1.
			cθ = 1
		end
		target_deffs[i] = rad2deg(acos(cθ))

		tar_dir_q = targetqua * cam_origindir / targetqua
		tar_dir[i, :] = [tar_dir_q.q1, tar_dir_q.q2, tar_dir_q.q3]

		# 撮影地点の画像上の軌跡
		if norm(i - shoot_time) < 120
			sat2tar_ecef = targetpos_ecef - (x_ecef_log[i, :] + x_ecef_error)
			sat2tar_seof = ecef_to_DCM(x_ecef_log[i, :] + x_ecef_error, v_ecef_log[i, :], true) * sat2tar_ecef
			sat2tar_seof = sat2tar_seof / norm(sat2tar_seof) # 衛星→撮影対象の単位方向ベクトル
			sat2tar_scsfqua = qua \ sat2tar_seof * qua
			sat2tar_scsf = [sat2tar_scsfqua.q1, sat2tar_scsfqua.q2, sat2tar_scsfqua.q3]
			sat2tar_scsf = sat2tar_scsf / norm(sat2tar_scsf)                  # カメラの単位方向ベクトル
			tarpos_log[j, :] = [sat2tar_scsf[1] * x_geod_log[i, 3] / sat2tar_scsf[3], sat2tar_scsf[2] * x_geod_log[i, 3] / sat2tar_scsf[3]]
			# if sat2tar_scsf[3] < 0
			# 	tarpos_log[j, :] = [NaN, NaN]
			# end
			

			#if norm(targetlocus_picture[j, 1]) > picture_xmax || norm(targetlocus_picture[j, 2]) > picture_ymax
			#	targetlocus_picture[j, : ] = [NaN, NaN]
			#end

			j = j + 1

			
			if norm(i-shoot_time)%12 == 0
				targetqua = targetqua_dicision(targetpos_ecef, x_ecef_log[i+12, :], v_ecef_log[i+12, :], sat_axisval)
			end
			
			# targetqua = targetqua_dicision(targetpos_ecef, x_ecef_log[i, :], v_ecef_log[i, :], sat_axisval)
		end

		



	end

	plot = true

	if plot == true
		plot_2scalar([1:DataNum],dotvs,"dotvs")
		plot_vec(direct_on_SCOFs,"direct_on_SCOFs")
		
		plot_2scalar([1:DataNum],atoms_denses,"atoms_dens")
		plot_2scalar([1:DataNum],x_geod_log[:,3],"height")
		plot_vec_range(mag_vecs,[-50000, 50000], [-50000, 50000], [-50000, 50000],"mag_vec")
		plot_3scalar([1:DataNum], mag_vecs[1:DataNum, 1], mag_vecs[1:DataNum, 2], mag_vecs[1:DataNum, 3], ["x", "y", "z"], "magvec_elements")
		plot_vec(sun_vecs,"sun_vec")
		plot_vec(x_ecef_log,"x_ecef_log")
		plot_vec(v_ecef_log,"v_ecef_log")
		plot_vec(x_geod_log,"x_geod_log")
		plot_2scalar(rad2deg.(x_geod_log[:,2]),rad2deg.(x_geod_log[:,1]),"x_geod_log_2d")
		plot_2scalar(rad2deg.(x_geod_log[shoot_time-120:shoot_time+120,2]),rad2deg.(x_geod_log[shoot_time-120:shoot_time+120,1]),"x_geod_log_2d_2")
		plot_vec(disturbance, "disturbance")
		plot_3scalar([1:DataNum], sat_ω[1:DataNum, 1], sat_ω[1:DataNum, 2], sat_ω[1:DataNum, 3], ["x", "y", "z"], "omega")
		plot_3scalar([1:DataNum], airtorques[1:DataNum, 1], suntorques[1:DataNum, 1], magtorques[1:DataNum, 1], ["air", "sun", "mag"], "x_torques")
		plot_3scalar([1:DataNum], airtorques[1:DataNum, 2], suntorques[1:DataNum, 2], magtorques[1:DataNum, 2], ["air", "sun", "mag"], "y_torques")
		plot_3scalar([1:DataNum], airtorques[1:DataNum, 3], suntorques[1:DataNum, 3], magtorques[1:DataNum, 3], ["air", "sun", "mag"], "z_torques")
		plot_3scalar([1:DataNum], T_reqs[1:DataNum, 1], T_reqs[1:DataNum, 2], T_reqs[1:DataNum, 3], ["x", "y", "z"], "request_torques")
		plot_vec(T_reqs, "request_torques_vec")
		plot_3scalar([1:DataNum], mtq_currentlog[:, 1], mtq_currentlog[:, 2], mtq_currentlog[:, 3], ["x", "y", "z"], "MTQcurrent")
		plot_3scalar([1:DataNum], sat_attqua_elements[1:DataNum, 2], sat_attqua_elements[1:DataNum, 3], sat_attqua_elements[1:DataNum, 4], ["q1", "q2", "q3"], "attqua_1")
		plot_2scalar([1:DataNum], sat_attqua_elements[1:DataNum, 1], "attqua_2")
		plot_vec(cam_dir, "cam_dir")
		plot_vecs(cam_dir, tar_dir, ["cam", "target"], "cam&tar_dir")
		plot_2scalar([1:DataNum], target_deffs, "target_deffs")
		plot_2scalar([shoot_time-120:shoot_time+120], target_deffs[shoot_time-120:shoot_time+120], "target_deffs2")
		# plot_2scalar(campos_log[:, 1], campos_log[:, 2], "campos")
		plot_2scalar(tarpos_log[:, 1], tarpos_log[:, 2], "tarpos")
		plot_2scalar([1:DataNum], controlmethod[:], "controlmethod")
		plot_3scalar([1:DataNum], Terrors[:, 1], T_reqs[1:DataNum, 1], magtorques[1:DataNum, 1], ["-", "reqs", "out"], "x_torques_check")
		plot_3scalar([1:DataNum], Terrors[:, 2], T_reqs[1:DataNum, 2], magtorques[1:DataNum, 2], ["-", "reqs", "out"], "y_torques_check")
		plot_3scalar([1:DataNum], Terrors[:, 3], T_reqs[1:DataNum, 3], magtorques[1:DataNum, 3], ["-", "reqs", "out"], "z_torques_check")

		# 撮影限界	画像サイズ4:3として範囲換算
		picture_radius = x_geod_log[shoot_time, 3] * tand(cam_viewangle/2)
		picture_xmax = 0.80 * picture_radius
		picture_ymax = 0.60 * picture_radius
		plot_2scalar_range(tarpos_log[:, 1], tarpos_log[:, 2], [-1*picture_xmax, picture_xmax], [-1*picture_ymax, picture_ymax], "targetlocus")
		plot_2scalar(tarpos_log[:, 1], tarpos_log[:, 2], "targetlocus2")

		# plot_2scalar(targetlocus_picture[:, 1], targetlocus_picture[:, 2], "targetlocus")
		
	end

	
end

main()