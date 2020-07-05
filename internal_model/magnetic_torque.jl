using LinearAlgebra  
"""
mag_tor = magneric_torque(i_m, B)

磁気トルカによるトルクを計算

# Arguments
- 'i_m'：各磁気トルカに流れる電流
- 'B'：磁束密度ベクトル＠SCSF

# Return
- 'mag_tor'：磁気トルカによる原点周りのトルク

"""
function magnetic_torque(i_m, B)
    μ = 5000 #透磁率
    n = 400 #巻き数
	L = 80  #コア長さ
	D = 10  #コア直径
	m_r = [0., 0., 0.]#[0.01, 0.01, 0.01] #残留磁気モーメント
	m = zeros(4,3)
	T = zeros(4,3)
	mag_tor = zeros(3)

	S = π*(D/1000)^2/4  # 断面積 [m]
	p = L/D
	μeff = 1 / (1/μ + (log(p)-1)/(p^2))  # コア実効比透磁率


	#磁気モーメント計算

    m[1,:] = μeff*n*S*i_m[1]*[1,0,0]
    m[2,:] = μeff*n*S*i_m[2]*[0,1,0]
	m[3,:] = μeff*n*S*i_m[3]*[0,0,1]
	m[4,:] = m_r

	#各軸周りのトルク計算
	for i=1:4
		T[i,:] = cross(m[i,:], B)
	end

	#全トルク計算
	for i=1:3
		for j=1:4
			global mag_tor[i] += T[j, i]
		end
	end

	return mag_tor
end