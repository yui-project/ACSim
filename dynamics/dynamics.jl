include("differencial.jl")
include("RungeKutta.jl")
"""
トルク`T`を与えられた剛体の回転を刻み時間`dt`秒ぶん計算する関数
# Arguments
- `q` ::Quaternion `t`秒時点の姿勢を表す単位クォータニオン 衛星軌道面座標系の基底をボディ座標系の基底に変換する
- `ω` ::Vector     `t`秒時点の角速度ベクトル　ボディ座標系基準
- `T` ::Vector     剛体に与えられるトルク　ボディ座標系基準
- `I` ::Matrix     慣性テンソル
- `dt`::           刻み時間 型指定してないけど良いですよね

# Returns
- `(qk/norm(qk),wk)` ::Tuple{Quaternion,Vector} `t+dk`秒時点の姿勢を表すクォータニオンと角速度ベクトルの組
"""
function dynamics(q::Quaternion,ω::Vector,T::Vector,I::Matrix,dt)
    qk=RK4((q,ω),dt)
    ωk=RK4((ω,T,I),dt) 
    return qk/norm(qk),ωk
end

