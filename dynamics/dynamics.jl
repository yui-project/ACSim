include("differencial.jl")
include("RungeKutta.jl")
"""
`dynamics(q::Quaternion,ω::Vector,T::Vector,I::Matrix,dt)`  
トルク`T`を与えられた剛体の回転を刻み時間`dt`秒ぶん計算する関数
## Arguments
- `q::Quaternion` :`t`秒時点の姿勢を表す単位クォータニオン 衛星軌道面座標系の基底をボディ座標系の基底に変換する
- `ω::Vector`     :`t`秒時点の角速度ベクトル　ボディ座標系基準
- `T::Vector`     :剛体に与えられるトルク　ボディ座標系基準
- `I::Matrix`     :慣性テンソル
- `dt::Number`          :刻み時間

## Returns
- `(qk/norm(qk), ωk)::Tuple{Quaternion,Vector}` :`(t+dt)`秒時点の姿勢を表す単位クォータニオンと角速度ベクトルの組

## Example
`qk,ωk=dynamics(q,ω,T,I,dt)`
"""
function dynamics(q::Quaternion,ω::Vector,T::Vector,I::Matrix,dt::Number)
    qk=RK4((q,ω),dt)
    ωk=RK4((ω,T,I),dt) 
    return qk/norm(qk),ωk
end

function dynamics(q::Quaternion,ω::Vector,T::Vector,I::Matrix,dt::Number)
    qk,ωk=RK4((q,ω),(T,I),dt) 
    return qk/norm(qk),ωk
end
