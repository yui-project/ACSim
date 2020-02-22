include("differencial.jl")
"""
`RK4(F,dt)`
Runge-Kutta法で刻み時間`dt`秒後のFの値を計算する関数
## Arguments
- `F::Any`: 計算したい時間の関数Fの`t`秒時点の値 関数`dif`で微分値を求められるもの
- `dt::`  :刻み時間
## Returns
- `dt`秒後のFの値
"""
function RK4(F,dt)
    f=dif(F)
    F=F[1]
    k1=dt*f
    k2=dt*(F+k1)
    k3=dt*(F+0.5(k1+k2))
    k4=dt*(F+k1+k3)
    return F+(k1+2k2+2k3+k4)/6
end
