using ReferenceFrameRotations
using LinearAlgebra


# x_ecef,v_ecefから軌道面座標系へのDCMを計算
function ecef_to_DCM(x_ecef,v_ecef,with_print=false)
    x_ecef_uniten = x_ecef / norm(x_ecef, 2)
    v_ecef_uniten = v_ecef / norm(v_ecef, 2)

    ez = - x_ecef_uniten
    ey = cross(v_ecef_uniten,x_ecef_uniten)
    ex = cross(ey,ez)

    C = inv_rotation(DCM([ex;ey;ez]))

    if with_print == true
        println("x_ecef_uniten:",x_ecef_uniten)
        println("v_ecef_uniten:",v_ecef_uniten)
        println("           ez:",ez)
        println("           ey:",ey)
        println("           ex:",ex)

        println("          DCM:",C)
    end
    return C
end

#=
# z軸、y軸からDCMを求める
function get_DCM_zy(z,y)
    ez = z / norm(z,2)
    ey = y / norm(y,2)

    ex = cross()
    ex = cross(ey,ez)

    C = inv_rotation(DCM([ex;ey;ez]))

    if with_print == true
        println("x_ecef_uniten:",x_ecef_uniten)
        println("v_ecef_uniten:",v_ecef_uniten)
        println("           ez:",ez)
        println("           ey:",ey)
        println("           ex:",ex)

        println("          DCM:",C)
    end
    return C
    
end
=#


#=
"""
# 方向余弦行列を用いてベクトルを変換
v = [2.0,0.0,0.0]
r = [1.0,1.0,1.0]

x = 1:10
@gp x x.^2 "w l tit 'Parabola'"
save("test2.gp")

C = getDCM(v,r)
println(C)
"""
=#