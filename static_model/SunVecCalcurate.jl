include("Coordinate_Transform.jl")
using SatelliteToolbox

"""
@fn     SunVecCal
@input  jd:計算したい時刻@ユリウス通日
@output 太陽方向ベクトル[m]@太陽中心黄道面基準慣性座標系
"""
function SunVecCal(jd::Number)
    
    return sun_position_i(jd)
end


julian_day = 2458749.815
println(SunVecCal(julian_day))