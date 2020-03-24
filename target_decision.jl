using SatelliteToolbox
using LinearAlgebra

"""
shootingtime = shootingtime_decision(satplace_plans, targetplace, currenttime, limit_time, dt)

撮影時刻を決定する関数

# Argments
 - `satplace_plans`：現在時刻~撮影期限までの衛星位置ベクトル@ECEF
 - `targetplace`：撮影対象の位置ベクトル@ECEF
 - `currenttime`：現在時刻（ループ回数）
 - `limit_time`：撮影期限までの時間 （現在時刻 + limit_time [s]までに撮影する）
 - `dt`：シミュレータループ時間間隔
# Return
- `shootingtime`：撮影時刻
- `shooting_vec`：撮影時衛星->撮影地点方向ベクトル@ECEF

# Conditions
 - `cam_viewangle`：カメラの視野角 [deg]
 - `sat_axisval`:衛星の撮影時姿勢の可動域 [deg]
"""

function shootingtime_decision(satplace_plans, targetplace, currenttime, limit_time, dt)
    cam_viewangle = 40
    sat_axisval = 10
    shootingtime = 0
    cθ_shooting = 0
    shootingtime = 0
    shooting_vec = zeros(1,3)

    cθ_lim = cosd(cam_viewangle/4 + sat_axisval)
    for i = 1:Int(limit_time/dt)
        vec_stot= satplace_plans(i, :) - targetplace
        cθ = vec_stot * satplace_plans(i, :)'
        cθ = cθ / (norm(vec_stot) * norm(satplace_plans(i, :)))

        if cθ < cθ_lim        # 撮影範囲外にいる時は除外
            continue
        elseif cθ > cθ_shooting
            cθ_shooting = cθ
            shootingtime = i
            shooting_vec = vec_stot
        end
    end

    shooting_vec = shooting_vec / norm(shooting_vec)

    return shootingtime, shooting_vec
end







    
end