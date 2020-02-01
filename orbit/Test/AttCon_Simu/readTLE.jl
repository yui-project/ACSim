#= 本プログラムファイルがあるフォルダ内にある、TLEを入力したtxtファイルを読み込み、
   TLEの要素ごとに分割するプログラム                                           =#

txtName="ISS_TLE.txt"
open("ISS_TLE.txt","r") do fp
    lines = readlines(fp)
    line1 = split(lines[2])
    line2 = split(lines[3])

    # 軌道計算に使う要素のみをGlobal変数として抽出    ##open-end 内で定義した変数はopen-endの外では使えないため
    global Epoch = parse(Float64, line1[4])
    print("Epoch:元期 : ")
    println(Epoch)
    global Incl = parse(Float64, line2[3])
    print("Incliation:軌道傾斜角 : ")
    println(Epoch)
    global RAAN = parse(Float64, line2[4])
    print("Right Ascension of Ascending Node:昇交点赤経 : ")
    println(RAAN)
    global Ecce = parse(Float64, line2[5])/(10^length(line2[5]))
    print("Eccentricity:離心率 : ")
    println(Ecce)
    global AP = parse(Float64, line2[6])
    print("Argument Perigee:近地点離隔 : ")
    println(AP)
    global MA = parse(Float64, line2[7])
    print("Mean Anomaly:平均近点角 : ")
    println(MA)
    # 通算周回数が増え10000回を超えると2行目の8要素と9要素が結合するため、例外処理
    if(size(line2)[1]<9)
        global MM = parse(Float64, line2[8][1:length(line2[8])-6])
        print("Mean Motion:平均運動 : ")
        println(MM)
        global RNE = parse(Int32, line2[8][length(line2[8])-5:length(line2[8])-1])
        print("Revolution Number at Epoch:通算周回数 : ")
        println(RNE)
    else
        global MM = parse(Float64, line2[8])
        print("Mean Motion:平均運動 : ")
        println(MM)
        global RevNum = parse(Int32, line2[9][1:length(line2[9])-1])
        print("Revolution Number at Epoch:通算周回数 : ")
        println(RNE)
    end
end

