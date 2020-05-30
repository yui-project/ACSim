using SatelliteToolbox
using Dates
using CSV
using DataFrames
using Plots
gr()

include("orbit/orbit.jl")
include("external_model/external_model.jl")

DataNum =1000 #シミュレータ反復回数
dt = 30 ##シミュレータの計算間隔 [s]
start_time = DateTime(2019, 12, 19, 3, 27, 10)	#シミュレート開始時刻
TLEFileName = "./orbit/ISS_TLE.txt"

# 進行方向とSCOFx軸とのずれ（内積）
dotvs = zeros(DataNum)

# 軌道計算については先に行い、全  時間分を配列に保存する
JD_log, x_ecef_log, v_ecef_log, x_geod_log, tles = orbit_cal(DataNum,dt,start_time,TLEFileName)

#println(JD_log)
#μ = 3.98600441* (10^14)
#a = (μ / ((tles[1].n * 2π / 86400)^2))^(1/3) 


#sun_angle_xp = satellite_sun_angle_earth_pointing(JD_log[1], a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [1, 0, 0])
#sun_angle_xm = satellite_sun_angle_earth_pointing(JD_log[1], a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [-1, 0, 0])
#sun_vec = zeros(length(sun_angle_xp), 3)

#for i = 1:length(sun_angle_xp)
#    if isnan(sun_angle_xp[i])
#        sun_vec[i, 1] = -1*cos(sun_angle_xm[i])
#    else
#        sun_vec[i, 1] = cos(sun_angle_xp[i])
#    end
#end

#sun_angle_yp = satellite_sun_angle_earth_pointing(JD_log[1], a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [0, 1, 0])
#sun_angle_ym = satellite_sun_angle_earth_pointing(JD_log[1], a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [0, -1, 0])
#for i = 1:length(sun_angle_yp)
#    if isnan(sun_angle_yp[i])
#        sun_vec[i, 2] = -1*cos(sun_angle_ym[i])
#    else
#        sun_vec[i, 2] = cos(sun_angle_yp[i])
#    end
#end

#sun_angle_zp = satellite_sun_angle_earth_pointing(JD_log[1], a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [0, 0, 1])
#sun_angle_zm = satellite_sun_angle_earth_pointing(JD_log[1], a, tles[1].e, tles[1].i, tles[1].M, tles[1].Ω, 1, [0, 0, -1])
#for i = 1:length(sun_angle_zp)
#    if isnan(sun_angle_zp[i])
#        sun_vec[i, 3] = -1*cos(sun_angle_zm[i])
#    else
#        sun_vec[i, 3] = cos(sun_angle_zp[i])
#    end
#end

sun_vec = sunvector_model(JD_log[1], tles)


df = DataFrame(A=sun_vec[:,2], B=sun_vec[:,3], C=sun_vec[:,4])
df = DataFrame(sun_vec)
df |> CSV.write("output.csv",delim=",", writeheader=false)


plot3d(sun_vec[:,2], sun_vec[:,3], sun_vec[:,4])
savefig("./figs/SunVecTest2.png")