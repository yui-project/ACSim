using SatelliteToolbox

"""
T = disturbance(sun_vecs, air_density, current_qua)

擾乱によるトルクの計算

# Arguments
- `sun_vecs`: 太陽方向ベクトル@SEOF
- `air_density`：大気密度
- `sat_velocity`:衛星速度ベクトル@seof
- `current_qua`: 現在の姿勢クォータニオン

# Returns
- `T`: 擾乱によるトルクの合計

"""
function disturbance(sun_vecs, air_density, sat_velocity, current_qua)
    sun_vecs = [0.05350930426409988, -0.8245359539307662, -0.5632736590971152]
    air_density = 1.8883960963389464e-12
    sat_velocity = [norm([2825.0313914335334, -4280.551397611556, -5681.579003761858]), 0., 0.]
    current_qua = SatelliteToolbox.Quaternion(cos(π/4), sin(π/4), 0, 0)
    #Ts = sun_pressure(sun_vecs, current_qua)
    #Ta = air_pressure(air_density, sat_velocity, current_qua)
    #println(Ts)
    #println(Ta)
    T = 0#Ts #+ Ta
    return T
end



"""
T = sun_pressure(sun_vecs, current_qua, sat_size, cm)

太陽輻射圧の計算

# Arguments
- `sun_vecs`: 太陽方向ベクトル@SEOF
- `current_qua`: 現在の姿勢クォータニオン
- `sat_size` : 衛星の形状（x, y, z）方向長さ[x,y,z][m]
- `cm` : 衛星体心から見た重心の位置ベクトル   [x,y,z][m]

# Returns
- `T`: 太陽輻射圧によるトルク
"""
function sun_pressure(sun_vecs, current_qua, sat_size, cm)
    Fe = 1353. #太陽定数 [W/m2]
    c = 299792458 # 光速 [m/s]
    ρs = 0.6 #鏡面反射割合
    ρd = 0.3 #散乱反射割合
    ρa = 0.1 #吸収割合


    #太陽方向ベクトルのBody座標系変換
    sunvecs_scsfqua = current_qua \ sun_vecs * current_qua
    #println(sunvecs_scsfqua)
    sunvecs_scsf = [sunvecs_scsfqua.q1, sunvecs_scsfqua.q2, sunvecs_scsfqua.q3]
    sunvecs_scsf_unit = sunvecs_scsf / norm(sunvecs_scsf)

    T = [0., 0., 0.]

    #x面に働く力
    if sunvecs_scsf[1] >= 0
        n = [1., 0., 0.]
        r = [sat_size[1]/2, 0., 0.] - cm
    else
        n = [-1., 0., 0.]
        r = [-sat_size[1]/2, 0., 0.] - cm
    end
    A = sat_size[2]*sat_size[3]
    P = Fe/c

    F = - P*A*dot(n,sunvecs_scsf_unit)* ((ρa+ρd)*sunvecs_scsf_unit + (2*dot(sunvecs_scsf_unit, n)ρs+2/3*ρd)*n)
    T = T + cross(r,F)

    # y面に働く力
    if sunvecs_scsf[2] >= 0
        n = [0., 1., 0.]
        r = [ 0., sat_size[2]/2, 0.] - cm
    else
        n = [0., -1., 0.]
        r = [0., -sat_size[2]/2, 0.] - cm
    end
    A = sat_size[1]*sat_size[3]
    
    F = - P*A*dot(n,sunvecs_scsf_unit)* ((ρa+ρd)*sunvecs_scsf_unit + (2*dot(sunvecs_scsf_unit, n)ρs+2/3*ρd)*n)
    T = T + cross(r,F)

    # z面に働く力
    if sunvecs_scsf[3] >= 0
        n = [0., 0., 1.]
        r = [0., 0., sat_size[3]/2] - cm
    else
        n = [0., 0., -1.]
        r = [0., 0., -sat_size[3]/2] - cm
    end
    A = sat_size[1]*sat_size[2]
    
    F = - P*A*dot(n,sunvecs_scsf_unit)* ((ρa+ρd)*sunvecs_scsf_unit + (2*dot(sunvecs_scsf_unit, n)ρs+2/3*ρd)*n)
    T = T + cross(r,F)

    return T

end




"""
T = air_pressure(density, current_qua)

空力トルクの計算

# Arguments
- `density`: 大気密度
- `vel_seof`：衛星の速度ベクトル＠SEOF
- `current_qua`: 現在の姿勢クォータニオン
- `sat_size` : 衛星の形状（x, y, z）方向長さ[x,y,z][m]
- `cm` : 衛星体心から見た重心の位置ベクトル   [x,y,z][m]
- `Cd` : 抗力係数（宇宙空間では通常2~3）

# Returns
- `T`: 空力トルク
"""
function air_pressure(density, vel_seof, current_qua, sat_size, cm, Cd)
    # 速度ベクトルのBody座標系への変換
    vel_scsfqua = current_qua \ vel_seof * current_qua
    vel_scsf = [vel_scsfqua.q1, vel_scsfqua.q2, vel_scsfqua.q3]


    T = [0., 0., 0.]

    #x面に働く力
    if vel_scsf[1] >= 0
        n = [1., 0., 0.]
        r = [sat_size[1]/2, 0., 0.] - cm
    else
        n = [-1., 0., 0.]
        r = [-sat_size[1]/2, 0., 0.] - cm
    end
    A = sat_size[2]*sat_size[3]
    

    F = -1/2 *Cd * density * dot(n, vel_scsf) * vel_scsf * A
    dT = cross(r, F)
    T = T + dT

    # y面に働く力
    if vel_scsf[2] >= 0
        n = [0., 1., 0.]
        r = [ 0., sat_size[2]/2, 0.] - cm
    else
        n = [0., -1., 0.]
        r = [0., -sat_size[2]/2, 0.] - cm
    end
    A = sat_size[1]*sat_size[3]

    F = -1/2 *Cd * density * dot(n, vel_scsf) * vel_scsf * A
    dT = cross(r, F)
    T = T + dT

    # z面に働く力
    if vel_scsf[3] >= 0
        n = [0., 0., 1.]
        r = [0., 0., sat_size[3]/2] - cm
    else
        n = [0., 0., 1.]
        r = [0., 0., -sat_size[1]/2] - cm
    end
    A = sat_size[1]*sat_size[2]
    
    F = -1/2 *Cd * density * dot(n, vel_scsf) * vel_scsf * A
    dT = cross(r, F)
    T = T + dT


    return T
    
end

## 以下，修正前の擾乱計算関数
#=
"""
T = sun_pressure(sun_vecs, current_qua)

太陽輻射圧の計算

# Arguments
- `sun_vecs`: 太陽方向ベクトル@SEOF
- `current_qua`: 現在の姿勢クォータニオン

# Returns
- `T`: 太陽輻射圧によるトルク
"""
function sun_pressure(sun_vecs, current_qua)
    sat_size = [0.1 , 0.1 , 0.1] #衛星各辺長さ [x,y,z][m]
    cm = [0.005, 0.005, 0.005] #衛星重心のずれ [x,y,z][m]
    I = 1353. #太陽定数 [W/m2]
    c = 299792458 # 光速 [m/s]
    ρs = 0.6 #鏡面反射割合
    ρd = 0.3 #散乱反射割合
    ρa = 0.1 #吸収割合


    #太陽方向ベクトルのBody座標系変換
    sunvecs_scsfqua = current_qua \ sun_vecs * current_qua
    #println(sunvecs_scsfqua)
    sunvecs_scsf = [sunvecs_scsfqua.q1, sunvecs_scsfqua.q2, sunvecs_scsfqua.q3]

    T = [0., 0., 0.]

    #x面に働く力
    if sunvecs_scsf[1] >= 0
        n = [1., 0., 0.]
        cf = [sat_size[1]/2, 0., 0.]
    else
        n = [-1., 0., 0.]
        cf = [-sat_size[1]/2, 0., 0.]
    end
    S = -sunvecs_scsf / norm(sunvecs_scsf)
    A = sat_size[2]*sat_size[3]
    P = I/c

    F =P*A*dot(n,S)*((ρa+ρd)*S+(2*ρs+2/3*ρd)*n)
    r = cf - cm

    T = T + cross(r,F)

    # y面に働く力
    if sunvecs_scsf[2] >= 0
        n = [0., 1., 0.]
        cf = [ 0., sat_size[2]/2, 0.]
    else
        n = [0., -1., 0.]
        cf = [0., -sat_size[2]/2, 0.]
    end
    S = -sunvecs_scsf / norm(sunvecs_scsf)
    A = sat_size[1]*sat_size[3]
    P = I/c

    F = P*A*dot(n,S)*((ρa+ρd)*S+(2*ρs+2/3*ρd)*n)
    r = cf - cm

    T = T + cross(r,F)

    # z面に働く力
    if sunvecs_scsf[3] >= 0
        n = [0., 0., 1.]
        cf = [0., 0., sat_size[3]/2]
    else
        n = [0., 0., - 1.]
        cf = [0., 0., -sat_size[3]/2]
    end
    S = -sunvecs_scsf / norm(sunvecs_scsf)
    A = sat_size[1]*sat_size[2]
    P = I/c

    F = P*A*dot(n,S)*((ρa+ρd)*S+(2*ρs+2/3*ρd)*n)
    r = cf - cm

    T = T + cross(r,F)

    return T

end

"""
T = air_pressure(density, current_qua)

空力トルクの計算

# Arguments
- `density`: 大気密度
- `vel_seof`：衛星の速度ベクトル＠SEOF
- `current_qua`: 現在の姿勢クォータニオン

# Returns
- `T`: 空力トルク
"""
function air_pressure(density, vel_seof, current_qua)
    sat_size = [0.1 , 0.1 , 0.1] #衛星各辺長さ [x,y,z][m]
    cm = [0.005, 0.005, 0.005] #衛星重心のずれ [x,y,z][m]
    Cd = 1.12 #抗力係数（各面を正方形と近似）

    # 速度ベクトルのBody座標系への変換
    vel_scsfqua = current_qua \ vel_seof * current_qua
    vel_scsf = [vel_scsfqua.q1, vel_scsfqua.q2, vel_scsfqua.q3]


    T = [0., 0., 0.]

    #x面に働く力
    if vel_scsf[1] >= 0
        n = [1., 0., 0.]
        cf = [sat_size[1]/2, 0., 0.]
    else
        n = [-1., 0., 0.]
        cf = [-sat_size[1]/2, 0., 0.]
    end
    A = sat_size[2]*sat_size[3]
    ξ = norm(cm)

    
    v_perp = dot(vel_scsf, -n) / norm(vel_scsf) * norm(vel_scsf)
    
    F = 1/2 * density * v_perp^2 * A * Cd * (-n)
    dT = cross(cf-cm, F)
    T = T + dT

    # y面に働く力
    if vel_scsf[2] >= 0
        n = [0., 1., 0.]
        cf = [ 0., sat_size[2]/2, 0.]
    else
        n = [0., -1., 0.]
        cf = [0., -sat_size[2]/2, 0.]
    end
    A = sat_size[1]*sat_size[3]

    v_perp = dot(vel_scsf, -n) / norm(vel_scsf) * norm(vel_scsf)
    
    F = 1/2 * density * v_perp^2 * A * Cd * (-n)
    dT = cross(cf-cm, F)
    T = T + dT

    # z面に働く力
    if vel_scsf[3] >= 0
        n = [0., 0., 1.]
        cf = [0., 0., sat_size[3]/2]
    else
        n = [0., 0., 1.]
        cf = [0., 0., -sat_size[1]/2]
    end
    A = sat_size[1]*sat_size[2]

    v_perp = dot(vel_scsf, -n) / norm(vel_scsf) * norm(vel_scsf)
    
    F = 1/2 * density * v_perp^2 * A * Cd * (-n)
    dT = cross(cf-cm, F)
    T = T + dT


    return T
    
end
=#