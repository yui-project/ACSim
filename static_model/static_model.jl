include("coordinate_transform.jl")
using SatelliteToolbox
using Dates

"""
static_model(datetime,r_ecef)

静的環境モデルに関わる全ての計算

# Arguments
- `datetime`: 時刻
- `r_ecef`: 衛星位置

# Returns
- `sun_vec`: 太陽方向ベクトル
- `shot_vec`: 撮影地点方向ベクトル
- `mag_vec`: 地磁場方向ベクトル
- `atoms_dens`: 大気密度スカラー

"""
function static_model(datetime,r_ecef)
    date = Date(datetime)
    jd = DatetoJD(date)

    mag_vel = mag_vec_cal(date,r_ecef)
    sun_vec = sun_vec_cal(jd)
    atoms_dens = atoms_dens_cal(height)
    shot_vec = shot_vec_cal(jd)
    return mag_vel, sun_vec, atoms_dens
end


"""
sun_vec_cal(jd::Number)

太陽方向ベクトルを算出する

# Arguments
- `jd`: 計算したい時刻@ユリウス通日

# Returns
- 太陽方向ベクトル[m]@太陽中心黄道面基準慣性座標系

# Examples
```jldoctest
julia> sun_vec_cal(2458749.815)
[-1.501323299678752e11, 3.644047456960307e6, 1.579694109989249e6]
```
"""
function sun_vec_cal(jd::Number)
    # 2019年秋分　julian_day = 2458749.815
    return sun_position_i(jd)
end

"""
@fn     shot_vec_cal
@input  撮影したい位置
@output 撮影地点方向ベクトル[m]@太陽中心黄道面基準慣性座標系
"""
function shot_vec_cal(jd::Number)
    # 2019年秋分　julian_day = 2458749.815
    return sun_position_i(jd)
end


"""
@fn     mag_vec_cal
@input  date:計算したい時刻
        r_ecef : 地心座標系(ECEF)による位置3行1列[m]
@output 
"""
function mag_vec_cal(date::Number, r_ecef::AbstractVector)
    # 筑波大 総合研究棟B r_ecef = [3.957729931663941e6,3.3091914769971482e6,3.737926200042794e6]
    r_geodetic = ECEFtoGeodetic(r_ecef)
    # r_geocentric = GeodetictoGeocentric(r_geocentric)
    return igrf12(date, norm(r_ecef), r_geodetic[1], r_geodetic[2], show_warns = true)
end


"""
@fn     atoms_dens_cal
@input  height:衛星高度[m]
@output density:待機密度
"""
function atoms_dens_cal(height)
    return expatmosphere(height)
end
