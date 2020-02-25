include("dynamics/dynamics.jl")
using Plots
#dynamics.jlのテスト
#Runge-Kutta法のテスト

len=63000
a=[zeros(len) zeros(len)]
a[1,2]=1.0
dt=0.01
#理論解
b=[sin.(range(0,(len-1)*dt,step=dt)) cos.(range(0,(len-1)*dt,step=dt))]

function dif(A::Array{Float64,1})
    θ=atan(A[1],A[2])
    cos(θ),-sin(θ)
end

for  i in 2:len
    a[i,:]=RK4((a[i-1,:]),dt)
    a[i,:]./=norm(a[i,:])
end