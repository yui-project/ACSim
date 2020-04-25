using SatelliteToolbox
using LinearAlgebra

#=
"""
shootingtime = shootingtime_decision(satplace_plans, targetplace, currenttime, limit_time, dt)

撮影時刻を決定する関数

# Argments
 - `satplace_plans`：現在時刻~撮影期限までの衛星位置ベクトル@ECEF
 - `targetplace`：撮影対象の位置ベクトル@ECEF
 - `currenttime`：現在時刻（ループ回数）
 - `limit_time`：撮影期限までの時間 （ループ回数　 現在時刻 ~ (currenttime + limit_time)*dt  の間で撮影する）
 - `dt`：シミュレータループ時間間隔
# Return
- `shootingtime`：撮影時刻
- `shooting_vec`：撮影時衛星->撮影地点方向ベクトル@ECEF

# Conditions
 - `cam_viewangle`：カメラの視野角 [deg]
 - `sat_axisval`:衛星の撮影時姿勢の可動域 [deg]
"""

function shootingtime_decision(satplace_plans, targetplace, currenttime, limit_time, cam_viewangle)
    cam_viewangle = 40
    sat_axisval = 30
    shootingtime = 0
    cθ_shooting = 0
    shootingtime = 0
    shooting_vec = zeros(1,3)

    cθ_lim = cosd(cam_viewangle/4 + sat_axisval)
    print("cθ_lim : ")
    println(cθ_lim)
    for i = 1:limit_time
        vec_stot= targetplace - satplace_plans[i, :]
        #=
        if i < 10
            println()
            print("vec_stot : ")
            println(vec_stot)
            print("targetplace : ")
            println(targetplace)
            print("satplace : ")
            println(satplace_plans[i , :])
        end
        =#

        cθ = dot(vec_stot, -satplace_plans[i, :])
        cθ = cθ / (norm(vec_stot) * norm(satplace_plans[i, :]))
        cθ_st = dot(targetplace, satplace_plans[i, :])
        cθ_st = cθ_st / (norm(targetplace) * norm(satplace_plans[i, :]))
        #=
        if i%10 == 0
            print("satplace : ")
            println(ECEFtoGeodetic(satplace_plans[i,:]))
            print("cθ : ")
            println(cθ)
            print("cθ_st : ")
            println(cθ_st)
        end
        =#
        if cθ < cθ_lim        # 撮影範囲外にいる時は除外
            continue
        elseif cθ_st < 0      # 地球の裏側にいる時は除外
            continue
        elseif cθ > cθ_shooting 
            cθ_shooting = cθ
            shootingtime = i
            shooting_vec = vec_stot
        end
    end

    shooting_vec = shooting_vec / norm(shooting_vec)

    print("cθ_shooting : ")
    println(cθ_shooting)

    return shootingtime, shooting_vec
end
=#

"""
shootingtime = shootingtime_decision2(satplace_plans, targetplace, currenttime, limit_time, cam_viewangle, sat_axisval)

撮影時刻を決定する関数

# Argments
 - `satplace_plans`：現在時刻~撮影期限までの衛星位置ベクトル@ECEF
 - `targetplace`：撮影対象の位置ベクトル@ECEF
 - `currenttime`：現在時刻（ループ回数）
 - `limit_time`：撮影期限までの時間 （ループ回数　 現在時刻 ~ (currenttime + limit_time)*dt  の間で撮影する）
 - `cam_viewangle`：カメラの視野角 [deg]
 - `sat_axisval`：衛星のポインティング時回転角度限界 [deg]
# Return
- `shootingtime`：撮影時刻(呼び出し元でのループ回数)
- `shooting_vec`：撮影時衛星->撮影地点方向ベクトル@ECEF
"""

function shootingtime_decision2(satplace_plans, targetplace, currenttime, limit_time,cam_viewangle, sat_axisval)
    shootingtime = 0

    shootingtime = 1
    satplace_geod = ECEFtoGeodetic(satplace_plans[1, :])
    targetplace_geod = ECEFtoGeodetic(targetplace)
    mindistance = sqrt((satplace_geod[1]-targetplace[1])^2+(satplace_geod[2]-targetplace[2])^2)
    
    # 衛星が撮影地点に最も近くなる時間を検索
    for i = 2:limit_time
        satplace_geod = ECEFtoGeodetic(satplace_plans[i, :])
        distance = sqrt((satplace_geod[1]-targetplace_geod[1])^2+(satplace_geod[2]-targetplace_geod[2])^2)
        if distance < mindistance
            shootingtime = i
            mindistance = distance
        end
    end

    #撮影範囲の限界値
    cθ_lim = cosd(cam_viewangle/4 + sat_axisval)
    
    # 撮影範囲内にいるか判断

    shooting_vec = targetplace - satplace_plans[shootingtime, :]
    shooting_vec = shooting_vec / norm(shooting_vec)
    print("shooting_vec : ")
    println(shooting_vec)
    print("satplace : ")
    println(satplace_plans[shootingtime, :])
    cθ = dot(shooting_vec, -satplace_plans[shootingtime, :])
    cθ = cθ / (norm(shooting_vec) * norm(satplace_plans[shootingtime, :]))
    
    if cθ < cθ_lim
        shooting_vec = [0., 0., 0.]
        shootingtime = 1
    else
        shootingtime = shootingtime + currenttime -1
    end

    r = sqrt(shooting_vec[1]^2+shooting_vec[2]^2)
    tφ = r/shooting_vec[3]
    if tφ > tand(sat_axisval)
        shooting_vec[1] = shooting_vec[1] * shooting_vec[3] * tand(sat_axisval) / r
        shooting_vec[2] = shooting_vec[2] * shooting_vec[3] * tand(sat_axisval) / r

        shooting_vec = shooting_vec / norm(shooting_vec)
    end
    
    return shootingtime, shooting_vec
end

"""
targetqua, rotqua = targetqua_dicision(shooting_vec, current_attqua)

撮影時の目標姿勢クォータニオンを求める

# Argments
 - `shooting_vec` : 撮影時の撮影地点方向単位ベクトル @SEOF
 - `current_attqua` : 現在の姿勢クォータニオン
 - `sat_axisval`：衛星のポインティング時回転角度限界 [deg]

# return
 - `targetqua` : 撮影時の目標姿勢クォータニオン
 - `rotqua` : 現在姿勢→目標姿勢になるために必要な回転を表すクォータニオン（SEOF基準）

# Conditions
 - `camera_initialdir` = 衛星が基準姿勢時の、カメラ撮影面の単位法線ベクトル @SEOF

"""
function targetqua_dicision(shooting_vec, current_attqua, sat_axisval)
    camera_initialdir = [0., 0., 1.]
    targetqua = SatelliteToolbox.Quaternion(0., 1., 0., 0.)
    
    c2θ = dot(camera_initialdir, shooting_vec) / (norm(camera_initialdir)*norm(shooting_vec))
    axis_seof = cross(camera_initialdir, shooting_vec)
    axis_seof = axis_seof / norm(axis_seof)

    cθ = sqrt((1+c2θ)/2)
    # θが回転限界を超えていたらθを回転限界まで抑える
    if cθ < cosd(sat_axisval)
        cθ = cosd(sat_axisval)
    end
    sθ = sqrt(1- cθ^2)

    targetqua = SatelliteToolbox.Quaternion(cθ, axis_seof[1]*sθ, axis_seof[2]*sθ, axis_seof[3]*sθ)
    rotqua = targetqua / current_attqua

    return targetqua, rotqua
end


"""
shootingplan_dicision
# 工事中
"""

function shootingplan_dicision(shootingtime, shooting_vec, current_time, current_attqua)
    targetqua = SatelliteToolbox.Quaternion(0., 1., 0., 0.)
    camera_initialdir = [0., 0., 1.]

    camera_dir = current_attqua * camera_initialdir / current_attqua

    θ = acos(dot(camera_dir, shooting_vec) / (norm(camera_dir)*norm(shooting_vec)))
    axis_seof = cross(camera_dir, shooting_vec)
    axis_seof = axis_seof / norm(axis_seof)

    axis_scsf = current_attqua \ axis_seof * current_attqua

end
