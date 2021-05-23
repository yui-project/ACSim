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
include("satellite/control_algorithm.jl")
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
	DataNum = 10000 #シミュレータ反復回数
	dt = 1 ##シミュレータの計算間隔 [s]
	start_time = DateTime(2019, 12, 19, 3, 27, 10)	#シミュレート開始時刻
	TLEFileName = "./orbit/ISS_TLE.txt"

	# 初期姿勢，角速度の設定
	sat_attqua_Initial = [cosd(15), sind(15)/sqrt(3), sind(15)/sqrt(3), sind(15)/sqrt(3)]
	sat_ω_Initial = [0., 0., 0.]

	# 撮影用パラメータの設定
	limit_time = DataNum                 # 計算開始時刻から "limit_time × dt" sec の間に撮影を行う
 	targetpos_geod = [0., -50., 25.7]    # 撮影対象の位置@Geodetic
	cam_viewangle = 40                   # カメラの視野角
	sat_axisval = 80                     # 撮影を許可する角度範囲（直下向きを0°とし，それと撮影時カメラ方向との角度差に対する制限）
	cam_origindir = [0., 0., 1.]         # カメラの方向ベクトル@SCSF
	picture_aspectratio = [4, 3]         # 撮影画像の縦横比
	num = 1                              # ターゲット軌跡記憶用配列のどこまでデータが入ったかを記憶する

	# 制御用パラメータの設定
	kp = 0.00000030                       # クロスプロダクト則比例ゲイン
	kr = 0.000030                         # クロスプロダクト則微分ゲイン
	mtq_maxcurrent = 0.002              # 磁気トルカの最大駆動電流
	mtq_scutter = 255                    # 磁気トルカの駆動電流分割数（" ± mtq_scutter" 段階で行う）
	Tmax = 1.0*10^(-7)                   # 出力トルクの最大値
	t_scatternum = 255                   # 出力トルクの分割数 (" ± t_scatternum" 段階で行う)
	I = [1/150 0.        0.;
		 0.        1/300 0.;
		 0.        0.        1/150]  # 衛星の慣性テンソル
	sat_size = [0.1 , 0.2 , 0.1]         # 衛星各辺長さ [x,y,z][m]
	cm = [0., 0., 0.]           # 衛星体心から見た重心の位置ベクトル [x,y,z][m]
	Cd = 2.0                             # 抗力係数（宇宙空間では通常2~3）
	target_updatefreq = 12               # 目標姿勢の更新頻度 [step/回]
	target_updaterange = 120             # 目標姿勢の更新を行う時間範囲（"撮影時刻 ± target_updaterange" の間は目標姿勢の更新を行う）
	CP2Bdot_delay = -1                    # CP制御の後 nステップはBdot制御を行わない
	CP2Bdot_count = 0                    # CP制御からのカウント数

	# 誤差検証用パラメータ
	x_ecef_error = [0., 0., 0. ]
	targetqua_error = SatelliteToolbox.Quaternion(cos(deg2rad(0)), 0, sin(deg2rad(0)), 0)


	#=
	データ保存用変数
	=#
	# 衛星位置，進行方向
	dotvs = zeros(DataNum)                    # 進行方向とSCOFx軸とのずれ（内積）
	direct_on_SCOFs = zeros(DataNum,3)        # SCOF上での衛星の進行方向
	# 外環境
	mag_vecs = zeros(DataNum,3)               # 地磁場@SEOF
	sun_vecs = zeros(DataNum,3)               # 太陽方向@ECEF
	atoms_denses = zeros(DataNum)             # 大気密度
	# 姿勢
	sat_attqua_elements = zeros(DataNum+1, 4) # 衛星姿勢クォータニオンの4要素
	sat_ω = zeros(DataNum+1, 3)               # 衛星の各軸回り角速度@SCSF
	sat_eularangle = zeros(DataNum, 3)      # 衛星姿勢（オイラー角）
	# トルク
	disturbance = zeros(DataNum, 3)           # 擾乱トルクの合計値@SCSF
	airtorques = zeros(DataNum, 3)            # 空力トルク

	
	suntorques = zeros(DataNum, 3)            # 太陽輻射圧トルク
	magtorques = zeros(DataNum, 3)            # 磁気トルク
	# 制御値
	mtq_currentlog = zeros(DataNum, 3)        # 磁気トルカ入力電流
	M_reqs = zeros(DataNum, 3)                # 所望磁気モーメント
	T_reqs = zeros(DataNum, 3)                # 所望トルク
	Terrors = zeros(DataNum, 3)               # 所望トルクの内，出せないトルク
	controlmethod = zeros(DataNum)            # どの制御を行っているかの判定
	# 撮影
	cam_dir = zeros(DataNum, 3)               # カメラ方向ベクトル@SEOF
	tar_dir = zeros(DataNum, 3)               # 目標方向ベクトル@SEOF
	target_diffs = zeros(DataNum)             # 目標方向ベクトルとカメラ方向ベクトルとの角度誤差
	tarpos_log = zeros(target_updaterange*2, 2)  # 撮影画像上でのターゲット位置の軌跡
	
	


	# 軌道・太陽方向計算については先に行い、全時間分を配列に保存する
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
		
		
	# 目標姿勢の計算
	targetpos_ecef = GeodetictoECEF(deg2rad(targetpos_geod[1]), deg2rad(targetpos_geod[2]), targetpos_geod[3])
	satqua = SatelliteToolbox.Quaternion(sat_attqua_elements[1,1], sat_attqua_elements[1,2], sat_attqua_elements[1,3], sat_attqua_elements[1,4])
	shoot_time = shootingtime_decision2(x_ecef_log[1:DataNum, :], targetpos_ecef, 1, limit_time, cam_viewangle, sat_axisval)
	println("shoottime", shoot_time)
	targetqua = targetqua_dicision(targetpos_ecef, x_ecef_log[shoot_time-target_updaterange, :], v_ecef_log[shoot_time-target_updaterange, :], sat_axisval)
	println("targetqua:", targetqua)
	targetqua = SatelliteToolbox.Quaternion(cosd(0), sind(0)/sqrt(3), sind(0)/sqrt(3), sind(0)/sqrt(3))
	
	
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
		mag_vec_seof = ecef_to_DCM(x_ecef_log[i,:],v_ecef_log[i,:],true)*mag_vec
		mag_vecs[i,:] = mag_vec_seof
		sun_vecs[i,:] = sun_vec
		atoms_denses[i] = atoms_dens
		println("         mag_vec:",mag_vec)
		println("         sun_vec:",sun_vec)
		println("      atoms_dens:",atoms_dens)
	

		direct_on_SCOFs[i,:] = ecef_to_DCM(x_ecef_log[i,:],v_ecef_log[i,:],true) * v_ecef_log[i,:]


		# 衛星姿勢
		qua = SatelliteToolbox.Quaternion(sat_attqua_elements[i,1], sat_attqua_elements[i,2], sat_attqua_elements[i,3], sat_attqua_elements[i,4])	

		
		#擾乱の計算
		v_scof = [norm(v_ecef_log[i,:]), 0., 0.]
		Ts = sun_pressure(sun_vecs[i,:], qua, sat_size, cm)
		Ta = air_pressure(atoms_denses[i], v_scof, qua, sat_size, cm, Cd)
		suntorques[i,:] = Ts
		airtorques[i,:] = Ta
		disturbance[i,:] = Ts + Ta


		# 地磁場ベクトル @SCSF の計算
		magvec_scsfqua = qua \ mag_vec_seof * qua
		magvec_scsf = [magvec_scsfqua.q1*10^(-9), magvec_scsfqua.q2*10^(-9), magvec_scsfqua.q3*10^(-9)]
		println("  scsf_mag_vec:", magvec_scsf)
		

		# 姿勢決定による衛星姿勢情報の取得
		qua_ad, vel_ad, omega_ad = attitude_empty(qua, v_ecef_log[i, :], sat_ω[i, :])

		#=
		# B-dot法による回転抑制制御
		if i==1
			M = B_dot(magvec_scsf, sat_ω[i, :], sat_ω[i, :])
		else
			M = B_dot(magvec_scsf, sat_ω[i, :], sat_ω[i-1, :])
		end
		=#

		
		
		# クロスプロダクト則による指向制御
		
		tarqua = targetqua_error * targetqua
		Treq, M = cross_product(tarqua,qua_ad, kp, kr, omega_ad, magvec_scsf)
		# Treq, M = pseudo_inverse(tarqua, qua_ad,kp, kr, omega_ad, magvec_scsf)
		# Treq, M, n = crossproduct_adj(targetqua, qua_ad, kp, kr, omega_ad, magvec_scsf)
		# controlmethod[i] = n

		#=
		# 所望トルクと地磁場のなす角が80°~100°の時はCrossProduct，それ以外はB-dot制御を行う
		tarqua = targetqua_error * targetqua
		Treq, M, n = cp_and_bdot(tarqua, qua_ad, omega_ad, kp, kr, magvec_scsf, boundangle)
		Terrors[i, :] = dot(Treq, magvec_scsf)/(norm(magvec_scsf)^2) * magvec_scsf
		=#
		
		
		M_reqs[i, :] = M
		println("request_Moment:", M)
		
		i_m = mm2current_theory(M)
		i_m = digitizing_3elesvec(i_m, mtq_maxcurrent, mtq_scutter, true)
		Tm = magnetic_torque(i_m, magvec_scsf)
		mtq_currentlog[i, :] = i_m
		magtorques[i,:] = Tm
		println("  current_mag:",i_m)
		
		# Treqの離散化
		Treq = digitizing_3elesvec(Treq, Tmax, t_scatternum, true)
		T_reqs[i, :] = Treq
		
		# Tm = collect(inv(ecef_to_DCM(x_ecef_log[i,:],v_ecef_log[i,:],true)) * Tm)
		
		# ダイナミクス
		# next_qua, ω = dynamics(qua, sat_ω[i,:], Ta+Ts+Tm, I, dt)
		# next_qua, ω = dynamics(qua, sat_ω[i,:], Tm, I, dt)
		# next_qua, ω = dynamics(qua, sat_ω[i,:], Treq, I, dt)
		next_qua, ω = dynamics(qua, sat_ω[i,:], Treq+Ta+Ts, I, dt)
		
		sat_attqua_elements[i+1, :] = [next_qua.q0, next_qua.q1, next_qua.q2, next_qua.q3]
		sat_ω[i+1, :] = ω
		println("quaternion:", qua)
		println("onega : ", ω)

		# カメラ方向
		cam_dir[i, :], target_diffs[i] = anglediffs(qua, targetqua, cam_origindir)

		tar_dir_q = targetqua * cam_origindir / targetqua
		tar_dir[i, :] = [tar_dir_q.q1, tar_dir_q.q2, tar_dir_q.q3]

		# オイラー角
		roll, pitch, yaw = Euler_Angles(qua \ targetqua)
		sat_eularangle[i, 1] = roll
		sat_eularangle[i, 2] = pitch
		sat_eularangle[i, 3] = yaw

		# 撮影地点の画像上の軌跡計算, 目標姿勢の更新
		if norm(i - shoot_time) < target_updaterange
			
			tarpos_log[num, :] = cal_targetlocus(targetpos_ecef, x_ecef_log, x_geod_log, v_ecef_log, x_ecef_error, qua, i)
			num = num + 1
			
			if norm(i-shoot_time) % target_updatefreq == 0
				# targetqua = targetqua_dicision(targetpos_ecef, x_ecef_log[i+target_updatefreq, :], v_ecef_log[i+target_updatefreq, :], sat_axisval)
			end
			
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
		plot_2scalar(rad2deg.(x_geod_log[shoot_time-target_updaterange:shoot_time+target_updaterange,2]),rad2deg.(x_geod_log[shoot_time-target_updaterange:shoot_time+target_updaterange,1]),"x_geod_log_2d_2")
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
		plot_vec(cam_dir[ Int(floor(DataNum*0.9)) : DataNum, : ], "cam_dir2")
		plot_vecs(cam_dir, tar_dir, ["cam", "target"], "cam&tar_dir")
		plot_2scalar([1:DataNum], target_diffs, "target_diffs")
		plot_2scalar([shoot_time-target_updaterange:shoot_time+target_updaterange], target_diffs[shoot_time-target_updaterange:shoot_time+target_updaterange], "target_diffs2")
		plot_2scalar(tarpos_log[:, 1], tarpos_log[:, 2], "tarpos")
		plot_2scalar([1:DataNum], controlmethod[:], "controlmethod")
		plot_3scalar([1:DataNum], Terrors[:, 1], T_reqs[1:DataNum, 1], magtorques[1:DataNum, 1], ["-", "reqs", "out"], "x_torques_check")
		plot_3scalar([1:DataNum], Terrors[:, 2], T_reqs[1:DataNum, 2], magtorques[1:DataNum, 2], ["-", "reqs", "out"], "y_torques_check")
		plot_3scalar([1:DataNum], Terrors[:, 3], T_reqs[1:DataNum, 3], magtorques[1:DataNum, 3], ["-", "reqs", "out"], "z_torques_check")
		
	

		plot_3scalar([1:DataNum], target_diffs, controlmethod*180, zeros(DataNum), ["targetdiffs", "controlmethod", "-"], "diffs_method_check")
		
		plot_3scalar([1:DataNum], sat_eularangle[:, 1], sat_eularangle[:, 2], sat_eularangle[:, 3], ["roll", "pitch", "yaw"], "EularAngleDiffs")
		
		# 撮影限界	画像サイズ4:3として範囲換算
		picture_radius = x_geod_log[shoot_time, 3] * tand(cam_viewangle/2)
		picture_xratio = picture_aspectratio[1] / norm(picture_aspectratio)
		picture_yratio = picture_aspectratio[2] / norm(picture_aspectratio)
		picture_xmax = picture_xratio * picture_radius
		picture_ymax = picture_yratio * picture_radius
		plot_2scalar_range(tarpos_log[:, 1], tarpos_log[:, 2], [-1*picture_xmax, picture_xmax], [-1*picture_ymax, picture_ymax], "targetlocus")
		

		# plot_2scalar(targetlocus_picture[:, 1], targetlocus_picture[:, 2], "targetlocus")
		
	end

	
end

main()