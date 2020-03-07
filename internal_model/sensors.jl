using LinearAlgebra

"""
test


sun_vol = sun_sensor(sun_vecs, sat_pos, sat_att)

太陽方向を電圧値として出力

# Argments
 - `sun_vecs`：太陽方向ベクトル
 - `ss_dir`：太陽センサ貼り付け面
 - `sat_att`：衛星の姿勢

# Return
- `sun_vol`：スポット光の位置(1*2)

# Conditions
 - `spot_height`：太陽センサのスポット - センサ間距離 [mm]
 - `sun_energy`：太陽光エネルギー密度@大気圏外 [kW/m2]
 - `ss_xLength`：太陽センサの光検出部長さ（x方向）[mm]
 - `ss_yLength`：太陽センサの光検出部長さ（y方向）[mm]
 - `ss_pointsensi`：太陽センサの光検出部位置分解能 [μm]
"""

function sun_sensor(sun_vecs, ss_dir)
	spot_height = 3.1
	ss_xLength = 9.0
	ss_yLength = 9.0
	ss_pointsensi = 1.5
	ssout_digits = 4
	sun_pos = zeros(1,2)
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
	# println(sunvecs_fix)
	# println(sunvecs_fix')

	if sunvecs_fix[3] >= 0
		# z成分 >= 0の時　センサに光が入射しないため出力なし
		sun_pos = [0.0 0.0]
	else
		# 太陽方向単位ベクトル、センサ-スポット間距離からセンサの出力値を計算
		sunvecs_norm = sunvecs_fix / norm(sunvecs_fix)
		sun_pos[1,1] = -1*sunvecs_norm[1]/sunvecs_norm[3]*spot_height
		sun_pos[1,2] = -1*sunvecs_norm[2]/sunvecs_norm[3]*spot_height
		# 光がセンサ検出外の時、例外処理
		if abs(sun_pos[1,1]) > ss_xLength || abs(sun_pos[1,2]) > ss_yLength
			sun_pos = [0.0 0.0]
		else
			sun_pos = round.(sun_pos ./ (ss_pointsensi/1000)) .* (ss_pointsensi/1000)
			sun_pos = round.(sun_pos, digits=ssout_digits)
		end
	end

	return sun_pos'
end

"""
mag_vol = mag_sensors(mag_vecs, sat_pos, sat_att)

地磁気ベクトルを電圧値として出力

# Argments
 - `mag_vecs`：地磁場ベクトル　[μT]
 - `output` : 出力方法（=`vol`：電圧, =`i2c`：I2Cシリアル通信）

# Return
 - `mag_vol`：地磁気ベクトルの電圧(1*3)
 - `mag_sel`：地磁気ベクトルのシリアル出力(1*3) 

# Conditions
 - `ms_range` : オフセットの調整なしで測定可能なセンサ動作範囲
 - `ms_sensi` : センサ分解能
 - `msout_minvol` : 出力最小値における電圧(output=volの時のみ)
 - `msout_maxvol` : 出力最大値における電圧(output=volの時のみ)
 - `msout_digit_vol`：電圧出力時の有効桁数
 - `msout_digit_i2c`：I2C出力時の有効桁数
 Conditions参考：BM1422AGMV
"""

function mag_sensor(mag_vecs, output)
	ms_range = 300   # ± μT
	ms_sensi = 0.042 # μT/LSB
	msout_minvol = 0.0 # V
	msout_maxvol = 5.0 # V
	msout_digit_vol = 6
	msout_digit_i2c = 3

	if output == "vol"
		msout_centvol = (msout_maxvol-msout_minvol)/2
		mag_vol = round.(mag_vecs ./ ms_sensi) .* (ms_sensi *(msout_centvol/2) / ms_range) .+ msout_centvol
		mag_vol = round.(mag_vol, digits=msout_digit_vol)
		return mag_vol

	elseif output == "i2c"
		mag_sel = round.(mag_vecs ./ ms_sensi) .* ms_sensi
		mag_sel = round.(mag_sel, digits=msout_digit_i2c)
		return mag_sel
	end

end