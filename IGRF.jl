using ReferenceFrameRotations
using LinearAlgebra
using Gnuplot


# g^m_n = gauss_coef_g[n][m+1]
gauss_coef_g = [-29442.0 -1501.0 00 ; -2445.1 3012.9 1676.7 0; 1350.7 -2352.3 1225.6 582.0]
# h^m_n = gauss_coef_h[n][m+1]
gauss_coef_h = [0.0 4797.1 0 0; 0.0 -2845.6 -641.9 0; 0.0 -115.3 244.9 -538.4]
# ルジェンドル
legendre = [0 0 0 0; 0 0 0 0; 0 0 0 0]

function cal_legendre(theta,N=3,M=4)
    # 初期値
    legendre_P[0][0] = 1
    legendre_P[1][0] = cos(theta) 

    for n in 2:N
        legendre_P[n][0] = ((2n -1) * cos(theta) * legendre_P[n-1][0] - (n-1) * legendre_P[n-2][0]) / n
        for m in 1:M


    return legendre
end

function cal_F()
    Xc = a


end


function cal_V(r,lambda,theta,a=6371.2)
    v = 0
    for n in 1:3
        for m in 1:4
            v += (a/r)^(n+1) * (gauss_coef_g[n][m] * cos(m * lambda) + gauss_coef_h[n][m] * sin(m * lambda)) * cal_legendre(cos(theta))
        end
    end
    v *= a
    return v
end

function cal_Xc()

    