using LinearAlgebra

"""
test

sun_vol = sun_sensor(sun_vecs, sat_pos, sat_att)

太陽方向を電圧値として出力

# Argments
 - `sun_vecs`：太陽方向ベクトル
 - `sat_pos`：衛星の位置
 - `sat_att`：衛星の姿勢

# Return
- `sun_vol`：太陽方向の電圧(1*2)
"""
function sun_sensor(sun_vecs, sat_pos, sat_att)
	sun_vol =  zeros(3)
	sun_vecs = zeros(3)
	sat_pos = zeros(3)
	sat_att = zeros(3)

	return sun_vol
end

"""
mag_vol = mag_sensors(mag_vecs, sat_pos, sat_att)

地磁気方向を電圧値として出力

# Argments
 - `mag_vecs`：地磁場方向ベクトル
 - `sat_pos`：衛星の位置
 - `sat_att`：衛星の姿勢

# Return
- `mag_vol`：地磁気方向の電圧(1*3)
"""
function mag_sensor(mag_vecs, sat_pos, sat_att)

	return mag_vol
end