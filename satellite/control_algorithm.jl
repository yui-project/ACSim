"""
digivec = digitizing_3elesvec(vec, max, steps, minus)

3要素ベクトルを離散化する．かつ，最大値を上回る場合は最大値に抑える

# Argument
 - `vec` : 離散化対象の3要素ベクトル
 - `max` : ベクトルの要素最大値
 - `steps` : 離散化段階（正負方向への離散化であれば，正方向に何段階分割するか） 
 - `minus` : 要素に負値を許容するか，許容するならtrue，しないなら非記入
# return
 - `digivec` : 離散化後のベクトル
"""

function digitizing_3elesvec(vec, max, steps, minus = false)
    resol = max / steps
    digivec = round.(vec / resol) * resol
    for i=1:3
        if norm(digivec[i]) > max
            digivec[i] = digivec[i] * max / norm(digivec[i])
        end

        if digivec[i] < 0. && minus == false
            digivec[i] = 0.
        end
    end

    return digivec
end


"""
Treq, M, n = cp_and_bdot(tarqua, qua_ad, omega_ad, kp, kr, magvec_scsf, boundangle)

クロスプロダクト則による所望トルクと地磁場のなす角が 90 ± n deg の時はクロスプロダクト則，他の時はBdot則による所望トルクに切り替えて制御を行う．

# Argument
 - `tarqua` : 目標姿勢クォータニオン
 - `qua_ad` : 現在姿勢クォータニオン
 - `omega_ad` : 現在の回転角速度
 - `kp` : クロスプロダクト則の比例ゲイン
 - `kr` : クロスプロダクト則の微分ゲイン
 - `magvec_scsf` : 地磁場ベクトル＠SCSF座標系
 - `boundangle` : クロスプロダクト則とB-dot則を切り替える角度（上記説明のnにあたる）
# Return
 - `Treq` : 所望トルク
 - `M` : 所望磁気モーメント
 - `n` : 制御方式（n=0でクロスプロダクト，n=1でBdot）

"""
function cp_and_bdot(tarqua, qua_ad, omega_ad, kp, kr, magvec_scsf, boundangle)
    n = 0
    Treq, M = cross_product(tarqua,qua_ad, kp, kr, omega_ad, magvec_scsf)
    φ = dot(Treq, magvec_scsf) / norm(Treq) / norm(magvec_scsf)
    if norm(φ) > cos(deg2rad(boundangle))
        if i==1
            M = B_dot(magvec_scsf, sat_ω[i, :], sat_ω[i, :])
        else
            M = B_dot(magvec_scsf, sat_ω[i, :], sat_ω[i-1, :])
        end

        Treq = cross(M, magvec_scsf)
        n = 1
        
    end
    
    return Treq, M, n
end

"""
vec1, θ = anglediffs(qua, targetqua, cam_origindir)


"""
function anglediffs(qua, targetqua, cam_origindir)
    qua1 = qua * cam_origindir / qua
    vec1 = [qua1.q1, qua1.q2, qua1.q3]
    qua2 = targetqua * cam_origindir / targetqua
    vec2 = [qua2.q1, qua2.q2, qua2.q3]
    cθ = dot(vec1, vec2) / norm(vec1) / norm(vec2)
    if cθ > 1.
        cθ = 1
    end
    θ = rad2deg(acos(cθ))
    
    return vec1, θ
end

"""
tarpos = cal_targetlocus(targetpos_ecef, x_ecef_log, x_geod_log, v_ecef_log, x_ecef_error, qua, i)
"""
function cal_targetlocus(targetpos_ecef, x_ecef_log, x_geod_log, v_ecef_log, x_ecef_error, qua, i)
    sat2tar_ecef = targetpos_ecef - (x_ecef_log[i, :] + x_ecef_error)
    sat2tar_seof = ecef_to_DCM(x_ecef_log[i, :] + x_ecef_error, v_ecef_log[i, :], true) * sat2tar_ecef
    sat2tar_seof = sat2tar_seof / norm(sat2tar_seof) # 衛星→撮影対象の単位方向ベクトル
    sat2tar_scsfqua = qua \ sat2tar_seof * qua
    sat2tar_scsf = [sat2tar_scsfqua.q1, sat2tar_scsfqua.q2, sat2tar_scsfqua.q3]
    sat2tar_scsf = sat2tar_scsf / norm(sat2tar_scsf)                  # カメラの単位方向ベクトル
    tarpos = [sat2tar_scsf[1] * x_geod_log[i, 3] / sat2tar_scsf[3], sat2tar_scsf[2] * x_geod_log[i, 3] / sat2tar_scsf[3]]

    return tarpos
end

    