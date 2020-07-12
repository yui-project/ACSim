include("../plot/plot_plots.jl")
# include("../plot/plot_makie.jl")

using LinearAlgebra


"""
matrix = cross_matrix(vector)

2つのベクトルの外積cross(a, v)のvを行列化し、v*aで外積を計算できるようにする

# Argments
 - `vector`：外積するベクトル

# Return
 - `matrix`：行列化したベクトル
"""
function cross_matrix(vector)
	matrix = zeros(3,3)
	matrix[1,1] = 0
	matrix[1,2] = vector[3]
	matrix[1,3] = -vector[2]
	matrix[2,1] = -vector[3]
	matrix[2,2] = 0
	matrix[2,3] = vector[1]
	matrix[3,1] = vector[2]
	matrix[3,2] = -vector[1]
	matrix[3,3] = 0
	return matrix
end


"""
t_req  = lyapunov_torque(sat_tar, sat_att, kp, kr, ω)

リアプノフ関数から目標姿勢までに必要なトルクを求める

# Argments
 - `sat_tar`：目標姿勢クォータニオン
 - `sat_att`：現在の衛星姿勢クォータニオン
 - `kp`：ポイントゲイン
 - `kr`：レートゲイン
 - `ω`：位置ベクトルの角速度（ベクトル）

# Return
 - `t_req`：必要トルク
"""
function lyapunov_torque(sat_tar, sat_att, kp, kr, ω)
	q_e = zeros(3)
	q = zeros(3)
	q_t = zeros(3)
	q0 = sat_att.q0
	q[1] = sat_att.q1
	q[2] = sat_att.q2
	q[3] = sat_att.q3
	q0_t = sat_tar.q0
	q_t[1] = sat_tar.q1
	q_t[2] = sat_tar.q2
	q_t[3] = sat_tar.q3
	#println("q", q)
	#println("qt", q_t)
	q_e = - q0*q + q0*q_t - cross(q,q_t)
	#println("qe", q_e)
	return kp*q_e - kr*ω
end


"""
m = B_dot()

B-dot法により
姿勢制御する際に必要な磁気モーメントを求める

# Argments
 - `B`：地磁気ベクトル@ECEF
 - `ω`：位置ベクトルの角速度（ベクトル)

# Return
- `m`：必要磁気モーメント
"""
function B_dot(B, ω, ω_b)
	k = zeros(3)

	k[1] = 10000
	k[2] = 10000
	k[3] = 10000

	m = -1 * k .* cross(B, ω)

	return m
end


"""
m = cross_product(sat_tar, sat_att, kp, kr, ω, B)

cross_product法により
姿勢制御する際に必要な磁気モーメントを求める

# Argments
 - `sat_tar`：目標姿勢クォータニオン
 - `sat_att`：現在の衛星姿勢クォータニオン
 - `kp`：ポイントゲイン
 - `kr`：レートゲイン
 - `ω`：位置ベクトルの角速度（ベクトル）
 - `B`：地磁気ベクトル@ECEF

# Return
 - `m`：必要磁気モーメント
"""
function cross_product(sat_tar, sat_att, kp, kr, ω, B)
	m = zeros(3)
	#必要トルクを求める
	t_req = lyapunov_torque(sat_tar, sat_att, kp, kr, ω)
	#println("t_req", t_req)
	#必要磁気モーメントを求める
	m = -(cross(t_req, B/norm(B)))/(norm(B))
	return t_req, m
end

"""
m = neo_cross_product(sat_tar, sat_att, kp, kr, ω, B)

cross_product法により
姿勢制御する際に必要な磁気モーメントを求める

# Argments
 - `sat_tar`：目標姿勢クォータニオン
 - `sat_att`：現在の衛星姿勢クォータニオン
 - `kp`：ポイントゲイン
 - `kr`：レートゲイン
 - `ω`：位置ベクトルの角速度（ベクトル）
 - `B`：地磁気ベクトル@ECEF

# Return
 - `m`：必要磁気モーメント
"""
function neo_cross_product(sat_tar, sat_att, kp, kr, ω, B)
	m = zeros(3)

	
	#必要トルクを求める
	t_req = lyapunov_torque(sat_tar, sat_att, kp, kr, ω)
	#println("t_req", t_req)
	#必要磁気モーメントを求める
	m = -(cross(t_req, B/norm(B)))/(norm(B))
	return t_req, m
end

"""
m = pseudo_inverse(sat_tar, sat_att, kp, kr, ω, B)

擬似逆行列により
姿勢制御する際に必要な磁気モーメントを求める

# Argments
 - `sat_tar`：目標姿勢クォータニオン
 - `sat_att`：現在の衛星姿勢クォータニオン
 - `kp`：ポイントゲイン
 - `kr`：レートゲイン
 - `ω`：位置ベクトルの角速度（ベクトル）
 - `B`：地磁気ベクトル@ECEF
 - `ω`：位置ベクトルの角速度

# Return
- `m`：必要磁気モーメント
"""
function pseudo_inverse(sat_tar, sat_att, kp, kr, ω, B)
	m = zeros(3,3)
	B_mat = zeros(3,3)
	B_mat = cross_matrix(B)
	#B_matの擬似逆行列
	B_pse = pinv(B_mat)
	#必要トルクを求める
	t_req = lyapunov_torque(sat_tar, sat_att, kp, kr, ω)
	m = B_pse*t_req
	return m
end


"""
out_cur = attitude_control()

目標姿勢まで姿勢制御する

# Argments
 - `sat_att`：現在の衛星姿勢クォータニオン
 - `sat_tar`：目標姿勢クォータニオン
 - `ω`：位置ベクトルの角速度（ベクトル）
 - `I`：慣性テンソル

# Return
- `out_cur`：加える電流値
"""
function attitude_control(sat_att, sat_tar, B, ω, I)
    
end


"""
i = mm2current_theory(M)

出力したい磁気モーメントから磁気トルカに流す電流値（理論値）を求める

#Argments
M : 必要とされる磁気モーメント [A/m2]

#Return
i : 磁気トルカに流す電流 [A]
"""
function mm2current_theory(M)
	n = 400  # 巻き数　[巻]
	L = 80  # コア長さ [mm]
	D = 10  # コア直径 [mm]
	μ = 5000  #コア初期比透磁率

	S = π*(D/1000)^2/4  # 断面積 [m]
	p = L/D
	μeff = 1 / (1/μ + (log(p)-1)/(p^2))  # コア実効比透磁率

	i = M / (μeff*n*S)

	return i
end


"""
i = mm2current_measure(M)

出力したい磁気モーメントから磁気トルカに流す電流値（実験値）を求める

#Argments
M : 必要とされる磁気モーメント [A/m2]

#Return
i : 磁気トルカに流す電流 [A]
"""
function mm2current_measure(M)
	slope = [0.92, 0.92, 0.92]
	intercept = [0.001, 0.001, 0.001]  # ix = slope[1] * Mx + intercept[1]となる様な切片、傾き

	i = slope .* M + intercept

	return i
end
