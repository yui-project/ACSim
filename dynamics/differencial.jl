using SatelliteToolbox
using LinearAlgebra

function dif((q,ω)::Tuple{Quaternion,Vector})
    0.5q*ω
end

function dif((ω,T,I)::Tuple{Vector,Vector,Matrix})
    (T-cross(ω,(I*ω)))*inv(I)
end