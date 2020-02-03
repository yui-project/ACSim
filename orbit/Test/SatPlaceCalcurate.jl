using SatelliteToolbox
#using Plots
#gr()

#=
設定パラメータ
=#
DataNum =100 #シミュレータ反復回数
dt = 5 ##シミュレータの計算間隔 [s]
t0 = DatetoJD(2019, 12, 19, 3, 27, 10)   #シミュレート開始時刻＠ユリウス日
TLEFileName = "./orbit/ISS_TLE.txt"

#事前準備関数定義
function set_SPCal(Name, t0)
    tles = read_tle(Name) #TLE読み込み
    #println("readTLE_Complete")
    orbp = init_orbit_propagator(Val{:sgp4}, tles[1]) #SGP4での軌道モデルの読み込み 
    #println("OrbitModelDownload_Complete")
    eop_IAU2000A = get_iers_eop(:IAU2000A)  #IAU-2000Aでの地球モデルの読み込み
    #println("EarthModelDownload_Complete")
    o, r, v = propagate_to_epoch!(orbp, t0 - dt/(24*60*60))
    #println("1stStepOrbitCalcurate_Complete")
    return orbp, eop_IAU2000A, o, r, v 
end

#ループ毎処理関数定義
function SPCal(n, t0, dt, orbp, eop_IAU2000A)
    #=  衛星位置の赤道面座標系での計算　=#
    o, r, v = step!(orbp, dt)
    #=  衛星位置の地心直交座標系への変換  =#
    xg = rECItoECEF(GCRF(), ITRF(), t0 + n*dt/(24*60*60), eop_IAU2000A) * r
    vg = rECItoECEF(GCRF(), ITRF(), t0 + n*dt/(24*60*60), eop_IAU2000A) * v
    return xg, vg
end


x_ecef_log = zeros(1000, 3)    #衛星位置＠地心直交座標系
x_geod_log = zeros(1000, 3)    #衛星位置＠測地座標系
v_ecef_log = zeros(1000, 3)    #衛星速度＠地心直交座標系


#シミュレーションループ
orbp, eop_IAU2000A, o, r, v = set_SPCal(TLEFileName, t0)


for i=1:DataNum
    x_sat, v_sat = SPCal(i, t0, dt, orbp, eop_IAU2000A)
    x_ecef_log[i,:] = x_sat
    v_ecef_log[i,:] = v_sat
    #println(x_ecef_log[i,:])
    print((i-1)*dt)
    print(", ")
    print(x_ecef_log[i,1])
    print(", ")
    print(x_ecef_log[i,2])
    print(", ")
    println(x_ecef_log[i,3])
    x_geo = ECEFtoGeodetic([x_sat[1]; x_sat[2]; x_sat[3]])
    x_geod_log[i,1]=x_geo[1]
    if x_geo[2]>0 
        x_geod_log[i,2]=x_geo[2]
    else
        x_geod_log[i,2]=x_geo[2]+2π
    end
    x_geod_log[i,2]=x_geo[2]
    x_geod_log[i,3]=x_geo[3]
    #println(x_geod_log[i,:])
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



#=***********************************************
緯度経度に変換してGifでグラフを作成したい場合はコレ
***********************************************=#
#=
anim = @animate for i=1:DataNum
    x_sat, v_sat = SPCal(i, t0, dt, orbp, eop_IAU2000A)
    x_ecef_log[i,:] = x_sat
    v_ecef_log[i,:] = v_sat
    x_geo = ECEFtoGeodetic([x_ecef_log[i,1]; x_ecef_log[i,2]; x_ecef_log[i,3]])
    x_geod_log[i,1]=x_geo[1]
    if x_geo[2]>0 
        x_geod_log[i,2]=x_geo[2]
    else
        x_geod_log[i,2]=x_geo[2]+2π
    end
    x_geod_log[i,3]=x_geo[3]
    println(i)
    plot(rad2deg.(x_geod_log[1:i,2]),rad2deg.(x_geod_log[1:i,1]))
end

gif(anim, "test5_anim.gif", fps = 10)
=#