using ReferenceFrameRotations
using LinearAlgebra
#using Gnuplot

# 位置と速度から座標変換のための方向余弦行列を求める
function get_DCM(v,r)
    e1 = v / norm(v, 2)
    e2 = - cross(r,v) / (norm(r,2) * norm(v,2))
    e3 = - r/ norm(r,2)
    C = inv_rotation(DCM([e1;e2;e3]))
    return C
end

"""
Satellite Centred Orbit Fixed
"""
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


function ecef_to_DCMa(x_ecef,v_ecef)
    C = get_DCM(v,r)
    return C * r
end



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