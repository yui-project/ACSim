using SatelliteToolbox
using Dates
using LinearAlgebra


"""
external_model(current_time,r_ecef, r_geod, cal_year=2020)

衛星外環境モデルに関わる全ての計算

# Arguments
- `current_time`: 現在時刻
- `r_ecef`: 衛星位置@ECEF
- `r_geod`: 衛星位置@Geodetic
- `cal_year`: 2020

# Returns
- `mag_vec_ecef`: 地磁場方向ベクトル@ECEF
- `sun_vecs`: 太陽方向ベクトル@SEOF
- `atoms_dens`: 大気密度スカラー

"""
function external_model(current_time,r_ecef, r_geod,r_eci_1st, eop_IAU2000A, sunvector_index, cal_year=2020)
    # 筑波大 総合研究棟B r_ecef = [3.957729931663941e6,3.3091914769971482e6,3.737926200042794e6]

    # ユリウス通日への変換
    jd = DatetoJD(current_time)

    # 地磁場ベクトル計算 @Geodetic to Geocentric
    # X軸：右手形をなすように設定
    # Y軸：東
    # Z軸：地球中心
    geoc_lat, geoc_h = GeodetictoGeocentric(r_geod[1], r_geod[3])
    mag_vec_geoc = igrf12(cal_year,geoc_h,geoc_lat,r_geod[2],Val{:geocentric})

    # mag_vec_geocのECEFへの変換
    # DCMの取得
    ez = - r_ecef/ norm(r_ecef,2) # 地球中心方向をz軸の単位ベクトルに
    ey = cross(ez,[0,0,1]) # ECEFのz軸との外積をy軸の単位ベクトルに
    ex = cross(ey,ez)
    C = DCM([ex;ey;ez])
    mag_vec_ecef = C * mag_vec_geoc


    # 太陽方向ベクトルの計算
     # 基準ベクトル(シミュレーション初期衛星位置@ECI)とのなす角を求める
    r_eci = rECEFtoECI(ITRF(), GCRF(), jd, eop_IAU2000A)*r_ecef
    cθ = dot(r_eci_1st,r_eci) / (norm(r_eci_1st)*norm(r_eci))
    sθ = norm(cross(r_eci_1st, r_eci)) / (norm(r_eci_1st)*norm(r_eci))
    if cθ > 1
        θ = π/2
    elseif cθ < -1
        θ = π/2
    else
        θ = acos(cθ)
    end
    if sθ < 0
        θ = 2π - θ
    end
     # indexからどの要素を抜き出すか
    index_num = round(Int64, 10*rad2deg(θ))
    sun_vecs = sunvector_index[index_num, 2:4]

    # 大気密度計算
    atoms_dens = expatmosphere(ECEFtoGeodetic(r_ecef)[3])

    return mag_vec_ecef, sun_vecs, atoms_dens
end

"""
sun_vec = sunvector_model(JD_log, tles)

太陽方向ベクトルの計算

# Arguments
- `JD_1st`: シミュレート開始時刻＠ユリウス日
- `tles`: シミュレートに用いるTLE情報 SatelliteToolbox.read_tle(Name)の出力形式

# Returns
- `sun_vec`: 衛星位置 - 太陽方向ベクトル@SEOF座標系 対応の行列
             1列目が衛星位置　（時刻JD_1stにおける衛星位置ベクトルとのなす角
             2-4列目がx,y,z方向成分

"""
function sunvector_model(JD_1st, tles)
    μ = 3.98600441* (10^14)  # 重力定数
    a = (μ / ((tles[1].n * 2π / 86400)^2))^(1/3)  # 軌道長半径 


     # x方向成分
    sun_angle_xp = satellite_sun_angle_earth_pointing(JD_1st, a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [1, 0, 0])
    sun_angle_xm = satellite_sun_angle_earth_pointing(JD_1st, a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [-1, 0, 0])
    sun_vec = zeros(length(sun_angle_xp), 4)

    #位置情報
    for i = 2:length(sun_angle_xp)
        sun_vec[i, 1] = sun_vec[i-1, 1] + 0.1
    end

    # x方向成分
    for i = 1:length(sun_angle_xp)
        if isnan(sun_angle_xp[i])
            sun_vec[i, 2] = -1*cos(sun_angle_xm[i])
        else
            sun_vec[i, 2] = cos(sun_angle_xp[i])
        end
    end

    # y方向成分
    sun_angle_yp = satellite_sun_angle_earth_pointing(JD_1st, a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [0, 1, 0])
    sun_angle_ym = satellite_sun_angle_earth_pointing(JD_1st, a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [0, -1, 0])
    for i = 1:length(sun_angle_yp)
        if isnan(sun_angle_yp[i])
            sun_vec[i, 3] = -1*cos(sun_angle_ym[i])
        else
            sun_vec[i, 3] = cos(sun_angle_yp[i])
        end
    end

    # z方向成分
    sun_angle_zp = satellite_sun_angle_earth_pointing(JD_1st, a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [0, 0, 1])
    sun_angle_zm = satellite_sun_angle_earth_pointing(JD_1st, a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [0, 0, -1])
    for i = 1:length(sun_angle_zp)
        if isnan(sun_angle_zp[i])
            sun_vec[i, 4] = -1*cos(sun_angle_zm[i])
        else
            sun_vec[i, 4] = cos(sun_angle_zp[i])
        end
    end

    return sun_vec

end
