#(n)STEP目の姿勢から(n+1)STEP目の姿勢を計算する
include("differencial.jl")
include("RungeKutta.jl")

function dynamics(q::Quaternion,ω::Vector,T::Vector,I::Matrix,dt)
    qk=RK4((q,ω),dt)
    ωk=RK4((ω,T,I),dt) 
    qk/norm(qk),ωk
end

