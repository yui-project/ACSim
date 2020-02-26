include("dynamics/dynamics.jl")
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

#クォータニオンとRK法による回転運動のシミュレーションのテスト(振り子)

l=10 #振り子の長さ
m=10 #重りの重さ
g=9.8*[0,-1,0] #重力加速度ベクトル
r0=l*[1,0,0] #振り子の初期位置ベクトル


#トルクを求める関数
function torque(q::Quaternion,(r0,m,g)::Tuple{Vector,Number,Vector})
    r=vect(q*r0*conj(q))
    T= cross(r,m*g)
    return T
end