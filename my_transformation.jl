using LinearAlgebra

"""
vout = eceftoseof(vin, satpos, satvel, vintype)

ECEFで表されるベクトルVinを、SEOFに座標変換する

# Argments
 - `vin`：座標変換を行う対象となるベクトル @ECEF
 - `satpos`：衛星位置ベクトル @ECEF
 - `satvel`：衛星速度ベクトル @ECEF
 - `vintype`：vinの種類（省略可、省略時vintype="dir"）　vintype="pos"->位置ベクトル vintype="dir"->方向ベクトル

# Return
 - `vout`：座標変換後のベクトル @SEOF

"""

function eceftoseof(vin, satpos, satvel, vintype="dir")
    # SEOFのz軸計算
    seof_z = -1.0 * satpos / norm(satpos)
    # SEOFのx軸計算
    normvel = satvel / norm(satvel)
    cθ = dot(seof_z, normvel)
    seof_x = normvel - cθ*seof_z
    # SEOFのy軸計算
    seof_y = cross(vec(seof_z), vec(seof_x))
    seof_y = seof_y / norm(seof_y)

    # ECEF -> SEOF()の座標変換行列の作成
    M_eceftoseof = seof_x'
    M_eceftoseof = hcat(M_eceftoseof, reshape(seof_y, 3, 1))
    M_eceftoseof = hcat(M_eceftoseof, seof_z')
    
    # 座標変換(回転)
    vout = M_eceftoseof \ vin'
    if vintype == "pos"
        satpos_seof = M_eceftoseof \ satpos'
        # 座標変換（原点の移動)
        vout = vout' .- satpos_seof'
    elseif vintype == "dir"
        vout = vout'
    else
        throw(ArgumentError("Enter \"dir\" or \"pos\" for vintype"))
    end

    return vout
    
end

"""
vout = eceftoseof(vin, satpos, satvel, vintype)

ECEFで表されるベクトルVinを、SEOFに座標変換する

# Argments
 - `vin`：座標変換を行う対象となるベクトル @ECEF
 - `satpos`：衛星位置ベクトル @ECEF
 - `satvel`：衛星速度ベクトル @ECEF
 - `vintype`：vinの種類（省略可、省略時vintype="dir"）　vintype="pos"->位置ベクトル vintype="dir"->方向ベクトル

# Return
 - `vout`：座標変換後のベクトル @SEOF

"""

function seoftoscsf(v, qua)
    vout = zeros(1,3)
    println(v)
    println(qua)
    vout_q = qua * vec(v) / qua

    vout[1]=vout_q.q1
    vout[2]=vout_q.q2
    vout[3]=vout_q.q3

    return vout
end

