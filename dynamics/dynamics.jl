#(n)STEP目の姿勢から(n+1)STEP目の姿勢を計算する
include("differencial.jl")
include("RungeKutta.jl")

function dynamics(q::Quaternion,ω::Vector,T::Vector,I::Matrix,dt)
    RK4((q,ω),dt),RK4((ω,T,I),dt)
end

