using LinearAlgebra
using SatelliteToolbox

"""
test


ss_out = sun_sensor(sun_vecs, ss_dir)

太陽方向を電流値として出力

# Argments
 - `sun_vecs`：太陽方向ベクトル
 - `ss_dir`：太陽センサ貼り付け面

# Return
- `ss_out`：出力端子から出力される電流値 [A]

# Conditions
 - `SUN_ENE_DENCE`：太陽光エネルギー密度 [W/m2]
 - `SPOT_HEIGHT`：太陽センサのスポット - センサ間距離 [mm]
 - `SPOT_DIA`：太陽センサのスポット径 [mm]
 - `sun_energy`：太陽光エネルギー密度@大気圏外 [kW/m2]
 - `SS_XLENGTH`：太陽センサの光検出部長さ（x方向）[mm]
 - `SS_YLENGTH`：太陽センサの光検出部長さ（y方向）[mm]
 - `SS_POINTSENSI`：太陽センサの光検出部位置分解能 [μm]
 - `SS_RESIST`：太陽センサの電極間抵抗 [Ω]
"""

function sun_sensor(sun_vecs, ss_dir)
	SUN_ENE_DENCE = 1400 
	SPOT_HEIGHT = 3.1
	SPOT_DIA = 1
	SS_XLENGTH = 9.0
	SS_YLENGTH = 9.0
	SS_POINTSENSI = 1.5
	SSOUT_DIGITS = 4
	SS_RESIST = 7000
	sun_pos = zeros(1,2)
	ss_out = zeros(1,4)
	sv_M = zeros(3,3)
	
	
	# 太陽センサ貼り付け面による座標変換
	if ss_dir == "x+"
		sv_M = [0.0 0.0 -1.0; 0.0 1.0 0.0; 1.0 0.0 0.0]
	elseif ss_dir == "x-"
		sv_M = [0.0 0.0 1.0; 0.0 1.0 0.0; -1.0 0.0 0.0]
	elseif ss_dir == "y+"
		sv_M = [1.0 0.0 0.0; 0.0 0.0 -1.0; 0.0 1.0 0.0]
	elseif ss_dir == "y-"
		sv_M = [1.0 0.0 0.0; 0.0 0.0 1.0; 0.0 -1.0 0.0]
	elseif ss_dir == "z+"
		sv_M = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
	elseif ss_dir == "z-"
		sv_M = [-1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 -1.0]
	else
		throw(ArgumentError(ss_dir, "Enter the surface attached a sun_sensor in ss_dir... x+,x-,y+,y-,z+, or z-"))
	end
	sunvecs_fix = sv_M*sun_vecs'

	if sunvecs_fix[3] >= 0
		# z成分 >= 0の時　センサに光が入射しないため出力なし
		sun_pos = [0.0 0.0]
	else
		# 太陽方向単位ベクトル、センサ-スポット間距離からセンサの出力値を計算
		sunvecs_norm = sunvecs_fix / norm(sunvecs_fix)
		sun_pos[1,1] = -1*sunvecs_norm[1]/sunvecs_norm[3]*SPOT_HEIGHT
		sun_pos[1,2] = -1*sunvecs_norm[2]/sunvecs_norm[3]*SPOT_HEIGHT
		# 光がセンサ検出外の時、例外処理
		if abs(sun_pos[1,1]) > SS_XLENGTH/2 || abs(sun_pos[1,2]) > SS_YLENGTH/2
			sun_pos = [0.0 0.0]
		else
			sun_pos = round.(sun_pos ./ (SS_POINTSENSI/1000)) .* (SS_POINTSENSI/1000)
			sun_pos = round.(sun_pos, digits=SSOUT_DIGITS)
			
			# 出力電流値の計算
			a = sqrt(SUN_ENE_DENCE * π * (SPOT_DIA/2/1000)^2 / SS_RESIST)
			b = 2*sun_pos[1,1] / SS_XLENGTH
			c = 2*sun_pos[1,2] / SS_YLENGTH
			if 1-c < 1-b
				ss_out[1] = a/2 - a*b/2 - a*c
			else
				ss_out[1] = a/2 - a*b - a*c/2
			end
			if ss_out[1] < 0
				ss_out[1] = 0.0
			end
			ss_out[2] = ss_out[1] + 1/2*a*(b+c)
			ss_out[3] = -1*ss_out[1] + 1/2*a - 1/2*a*b
			ss_out[4] = -1*ss_out[1] + 1/2*a - 1/2*a*c		
		end
	end

	return ss_out
end

"""
mag_vol = mag_sensor(mag_vecs)

地磁気ベクトルをシリアル出力

# Argments
 - `mag_vecs`：地磁場ベクトル　[μT]

# Return
 - `ms_out`：地磁気ベクトルのシリアル出力(1*3) 

# Conditions
 - `MS_RANGE` : オフセットの調整なしで測定可能なセンサ動作範囲
 - `MS_SENSI` : センサ分解能
 - `MSOUT_DIGIT`：I2C出力時の有効桁数
 Conditions参考：BM1422AGMV
"""

function mag_sensor(mag_vecs)
	MS_RANGE = 300   # ± μT
	MS_SENSI = 0.042 # μT/LSB
	MSOUT_DIGIT = 3

	ms_out = round.(mag_vecs ./ MS_SENSI) .* MS_SENSI
	ms_out = round.(ms_out, digits=MSOUT_DIGIT)

	return ms_out

end


"""
acce_out, gyro_out = gyro_sensor(v_ecef_pu, attqua1, attqua2, dt)

角速度、加速度ベクトルをシリアル出力

# Argments
 - `v_ecef_pu`：ループi-1,i回目の衛星位置@ECEFログ(3*2)
 - `attqua1`：ループi-2回目の衛星姿勢クォータニオン
 - `attqua1`：ループi-1回目の衛星姿勢クォータニオン
 - `dt`：シミュレータ計算時間間隔

# Return
 - `acce_out`：加速度ベクトル@ECEF
 - `gyro_out`：角速度ベクトル@Body

# Conditions
 - `GS_ACCERANGE` : 加速度センサ測定範囲 ±g
 - `GS_GYRORANGE` : 角速度センサ測定範囲　deg/s
 - `GS_ACCEBITS` : 加速度センサ1軸あたり使用ビット数
 - `GS_GYROBITS` : 角速度センサ1軸あたり使用ビット数
 - `GS_ACCEDIGITS` : 加速度センサ出力値有効桁数 点以下n桁
 - `GS_GYRODIGITS` : 角速度センサ出力値有効桁数 点以下n桁
 Conditions参考：BMI160
"""

function gyro_sensor(v_ecef_pu, attqua1, attqua2, dt)
	GS_ACCERANGE = 4
	GS_GYRORANGE = 250
	GS_ACCEBITS = 65536
	GS_GYROBITS = 65536
	GS_ACCEDIGITS = 4
	GS_GYRODIGITS = 4

	accerange = GS_ACCERANGE * 9.81
	accesensi = accerange*2 / GS_ACCEBITS

	gyrorange = GS_GYRORANGE / 180 * π
	gyrosensi = gyrorange*2 / GS_GYROBITS
	
	acce = (v_ecef_pu[2][:] - v_ecef_pu[1][:])/dt
	gyro = gyro_calcurate(attqua1, attqua2, dt)

	acce_out = round.(acce ./ accesensi) .* accesensi
	acce_out = round.(acce_out, digits=GS_ACCEDIGITS)
	gyro_out = round.(gyro ./ gyrosensi) .* gyrosensi
	gyro_out = round.(gyro_out, digits=GS_GYRODIGITS)
	return acce_out, gyro_out
	
end


"""
gyro = gyro_calcurate(qua1, qua2, dt)

状態1から状態2に一定速度で回転して変化する時の角速度ベクトルを求める
但し、dtのうちに ±180deg以上回転する場合は±180degに収まる値に変換される（+200deg -> -160deg, +380deg -> +20deg）
(±180deg以上/以下の判定が出来ていないため)

# Argments
 - `qua1`：状態1（変化前）の衛星姿勢を示すクォータニオン
 - `qua2`：状態2（変化後）の衛星姿勢を示すクォータニオン
 - `dt`：状態1 -> 状態2 に遷移するのに要する時間

# Return
 - `gyro`：角速度ベクトル＠Body座標系

# Conditions
"""

function gyro_calcurate(qua1, qua2, dt)
	ω=zeros(3,1)

	# qua1 -> qua2 へ変化するためのクォータニオン@ECEF (cq = Change Quaternion)を求める
	cq_ecef = qua2 / qua1
	if cq_ecef.q0 < 0
		cq_ecef = -1*cq_ecef
	end	

	# cq_ecefをオイラー軸/オイラー角に分離、オイラー軸をBody座標系に変換してから再結合、ついでにcqの時間微分（オイラー角/dt）を求める
	θ = acos(cq_ecef.q0)
	if sin(θ) != 0.0  # 例外処理：sinθ=0であればqua1=qua2 => 回転していない =>　ω=[0., 0., 0.]
		n = [cq_ecef.q1 / sin(θ); cq_ecef.q2 / sin(θ); cq_ecef.q3 / sin(θ)]
		dθ = θ/dt
		axe = qua2 * n / qua2
		dq = SatelliteToolbox.Quaternion(cos(dθ), axe.q1*sin(dθ), axe.q2*sin(dθ), axe.q3*sin(dθ))
		cq_body = SatelliteToolbox.Quaternion(cos(θ), axe.q1*sin(θ), axe.q2*sin(θ), axe.q3*sin(θ))
		
		# dq と角速度ベクトルの関係式より、角速度ベクトルを求める
		omega = 2* cq_body \ dq

		ω[1] = omega.q1
		ω[2] = omega.q2
		ω[3] = omega.q3
	end

	return ω
end
