include("dynamics.jl")
using Plots
using ReferenceFrameRotations
#dynamics.jlのテスト
#Runge-Kutta法のテスト

len=63000
a=[zeros(len) zeros(len)]
a[1,2]=1.0
dt=0.01
#理論解
b=[sin.(range(0,(len-1)*dt,step=dt)) cos.(range(0,(len-1)*dt,step=dt))]
#微分する関数
function dif(v::Array{Float64,1},A)
    θ=atan(v[1],v[2])
    cos(θ),-sin(θ)
end
#積分
for  i in 2:len
    a[i,:]=RK4((a[i-1,:]),0,dt)
    a[i,:]./=norm(a[i,:])
end