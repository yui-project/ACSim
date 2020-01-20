#include("orbit/orbit.jl")
include("static_model/static_model.jl")
#include("dynamic_model/dynamic_model.jl")
#include("dynamics/dynamics.jl")
include("satellite/satellite.jl")
using Dates


"""
main

メイン、以下嘘800

# Arguments
- `datetime`: 時刻
- `r_ecef`: 衛星位置

# Returns
- `sun_vec`: 太陽方向ベクトル
- `shot_vec`: 撮影地点方向ベクトル
- `mag_vec`: 地磁場方向ベクトル
- `atoms_dens`: 大気密度スカラー

"""
function main()
    
end