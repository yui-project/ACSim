using SatelliteToolbox
using LinearAlgebra

"""
`dif((q,ω)::Tuple{Quaternion,Vector})`
クォータニオンを微分する関数
## Argument
- `(q,ω)::Tuple{Quaternion,Vector}` :時刻`t`における姿勢を表すクォータニオンと角速度ベクトルの組

## Returns
- 時刻`t`における姿勢を表すクォータニオンを時間微分したもの
"""
function dif((q,ω)::Tuple{Quaternion,Vector})
    0.5q*ω
end

"""
`(ω,T,I)::Tuple{Vector,Vector,Matrix}`
角速度ベクトルを微分する関数
## Argument
- `(ω,T,I)::Tuple{Vector,Vector,Matrix}` :時刻`t`における角速度ベクトルと、剛体に加わるトルク、剛体の慣性テンソルの組

## Returns
- 時刻`t`における角速度ベクトルを時間微分したもの
"""
function dif((ω,T,I)::Tuple{Vector,Vector,Matrix})
    (T-cross(ω,(I*ω)))*inv(I)
end