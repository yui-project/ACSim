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
    μ = 1 #magnetic permeability
    n = 1 #the number of turns
    S = 1 #section area
	m_r = [1,1,1] #残留磁気モーメント
	m = zeros(4)
	T = zeros(4,3)
	mag_tor = zeros(3)

	#磁気モーメント計算
    m[1] = μ*n*S*i_m[1]*[1,0,0]
    m[2] = μ*n*S*i_m[2]*[0,1,0]
	m[3] = μ*n*S*i_m[3]*[0,0,1]
	m[4] = m_r

	#各軸周りのトルク計算
	for i=1:4
		T[i,:] = cross(m[i], B)
	end

	#全トルク計算
	for i=1:3
		for j=1:4
			global mag_tor[i] += T[j, i]
		end
	end

	return mag_tor
end