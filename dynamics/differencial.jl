using SatelliteToolbox
using LinearAlgebra

"""
`dif((q,ω)::Tuple{Quaternion,Vector})`
クォータニオンと角速度ベクトルを微分する関数
## Argument
- `(q,ω)::Tuple{Quaternion,Vector}` :時刻`t`における姿勢を表すクォータニオンと角速度ベクトルの組
- `(T,I)::Tuple{Vector,Matrix}` :時刻`t`におけるトルクベクトルと慣性テンソルの組
## Returns
- `(dq,dω)::Tuple{Quaternion,Vector}`時刻`t`における姿勢を表すクォータニオンと角速度ベクトルを時間微分したものの組
"""
function dif((q,ω)::Tuple{Quaternion,Vector},(T,I)::Tuple{Vector,Matrix})
    dq = 0.5q*ω
    dω = inv(I)*(T-cross(ω,(I*ω)))
    return (dq,dω)
end