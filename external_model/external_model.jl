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
- `sun_vec`: 太陽方向ベクトル@
- `mag_vec`: 地磁場方向ベクトル@Geodetic
- `atoms_dens`: 大気密度スカラー@

"""
function external_model(current_time,r_ecef, r_geod, cal_year=2020)
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



    # 太陽方向計算
    sun_vec = sun_position_i(jd)
    #satellite_sun_angle_earth_pointing(jd,)

    # 大気密度計算
    atoms_dens = expatmosphere(ECEFtoGeodetic(r_ecef)[3])

    return mag_vec_ecef, sun_vec, atoms_dens
end