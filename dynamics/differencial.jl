using SatelliteToolbox
using LinearAlgebra

function dif((q::Quaternion,ω::Vector))
    0.5q*ω
end

function dif((ω::Vector,T::Vector,I::Matrix))
    (T-cross(ω,(I*ω)))*inv(I)
end