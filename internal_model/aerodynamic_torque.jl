using LinearAlgebra
"""
aer_tor = aerodynamics_torque(C_d, rho, v, a, e_r, e_a, p_c)

空力トルクの計算

# Arguments
- `C_d`：抵抗係数
- `rho`：大気密度
- `v`：大気に対する衛星の相対速度
- `a`：衛星の各面の面積
- `e_r`：衛星の進行方向の単位ベクトル@SCSF（body）
- `e_a`：各面の法線ベクトル@SCSF（body）
- `p_c`：圧力中心の位置ベクトル@SCSF（body）

# Return
- `aer_tor`：空力トルク（1*3）

"""
function aerodynamic_torque(sur_num, C_d, rho, v, a, e_r, e_a, p_c)
	tor = zeros(3)
	cos = zeros(6)
	for i=1:sur_num
		#e_rとe_aのなす角の余弦を計算
		#cos(theta)>0の面が空気抵抗をうける面、cos(theta)<0の面が空気抵抗を受けない面
		ea = e_a[i,:]
		cos[i] = dot(e_r, ea)
		if cos[i]>0.0
			pc = p_c[i,:]
			tor += a[i]*cos[i]*cross(pc, e_r)
		end
	end
	aer_tor = -(1/2)*C_d*rho*v*v*tor

	return aer_tor
end