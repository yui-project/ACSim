include("differencial.jl")
"""
`RK4(F,A,dt)`
Runge-Kutta法で刻み時間`dt`秒後のFの値を計算する関数
## Arguments
- `F::Any`: 計算したい時間`t`の関数の`t`秒時点の値 関数`dif`で微分値を求められるもの
- `A::Any`:微分計算に必要な定数の組
- `dt::Number`  :刻み時間
## Returns
- `dt`秒後のFの値
"""
function RK4(F,A,dt)
    f=dif(F,A)
    k1=dt.*f
    k2=dt.*dif((F.+k1),A)
    k3=dt.*dif((F.+0.5.*(k1.+k2)),A)
    k4=dt.*dif((F.+k1.+k3),A)
    return F.+(k1.+2.0.*k2.+2.0.*k3.+k4)./6
end
