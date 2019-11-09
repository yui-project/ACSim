using SatelliteToolbox
using Plots
gr()

#=
設定パラメータ
=#
DataNum =1  #シミュレータ反復回数
dt = 60*5 ##シミュレータの計算間隔 [s]
t = DatetoJD(2019, 11, 07, 16, 10, 00)   #シミュレート開始時刻＠ユリウス日


#開始時のみ処理
tles = read_tle("ISS_TLE.txt") #TLE読み込み
orbp = init_orbit_propagator(Val{:sgp4}, tles[1]) #SGP4での軌道モデルの読み込み 
eop_IAU2000A = get_iers_eop(:IAU2000A)  #IAU-2000Aでの地球モデルの読み込み
epoch_JD = DatetoJD(2000+tles[1].epoch_year, 1, 1, 00, 00, 00) + tles[1].epoch_day  #TLEの元期＠ユリウス日 
o, r, v = propagate!(orbp, t - epoch_JD)

#ループ毎処理関数定義
function SatPlaceCal(n)
    #=  衛星位置の赤道面座標系での計算　=#
    o, r, v = step!(orbp, dt)
    #=  衛星位置の地心直交座標系への変換  =#
    xg = rECItoECEF(GCRF(), ITRF(), t + n*dt/(24*60*60), eop_IAU2000A) * r
    vg = rECItoECEF(GCRF(), ITRF(), t + n*dt/(24*60*60), eop_IAU2000A) * v
    return xg, vg
end

x_ecef_log = zeros(100, 3)    #衛星位置＠地心直交座標系
x_geod_log = zeros(100, 3)    #衛星位置＠測地座標系
v_ecef_log = zeros(100, 3)    #衛星速度＠地心直交座標系

#シミュレーションループ
for i=1:DataNum
    x_sat, v_sat = SatPlaceCal(i)
    x_ecef_log[i,:] = x_sat
    v_ecef_log[i,:] = v_sat
    x_geo = ECEFtoGeodetic([x_sat[1]; x_sat[2]; x_sat[3]])
    x_geod_log[i,1]=x_geo[1]
    x_geod_log[i,2]=x_geo[2]
    x_geod_log[i,3]=x_geo[3]
end



#=
#=************************************************
地心直交座標系での衛星位置を3Dプロットしたい場合はコレ
************************************************=#
println(x_ecef_log)

plot(x_ecef_log[:,1], x_ecef_log[:, 2], x_ecef_log[:, 3])
xlabel!("x")
ylabel!("y")
#zlabel!("z")

savefig("test2.png")
=#

#=
#=***********************************************
緯度経度に変換してGifでグラフを作成したい場合はコレ
***********************************************=#
anim = @animate for i=1:1
    x_sat, v_sat = SatPlaceCal(i)
    x_ecef_log[i,:] = x_sat
    v_ecef_log[i,:] = v_sat
    x_geo = ECEFtoGeodetic([x_ecef_log[i,1]; x_ecef_log[i,2]; x_ecef_log[i,3]])
    x_geod_log[i,1]=x_geo[1]
    x_geod_log[i,2]=x_geo[2]
    x_geod_log[i,3]=x_geo[3]
    plot(2*x_geod_log[1:i,2],x_geod_log[1:i,2])
end

gif(anim, "test4_anim.gif", fps = 10)
=#