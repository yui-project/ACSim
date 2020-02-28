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
d=1 #重りの直径

g=9.8*[0,0,-1] #重力加速度ベクトル
r0=l*[1,0,0] #振り子の初期位置ベクトル
#慣性テンソル
Ix=0.4m*(d^2)
Iy=0.4m*(d^2)+m*(l^2)
Iz=0.4m*(d^2)+m*(l^2)
I=
[Ix 00 00
;00 Iy 00
;00 00 Iz]


#トルクを求める関数
function torque(q::Quaternion,(r0,m,g)::Tuple{Vector,Number,Vector})
    r=vect(q*r0*conj(q))
    T= cross(r,m*g)
    return T
end

#力学的エネルギーを求める関数
function energy(q::Quaternion,ω::Vector,(r0,m,g)::Tuple{Vector,Number,Vector})
    r=vect(q*r0*conj(q))
    v=cross(r,ω)
    f=m*g
    E=0.5m*dot(v,v)-dot(r,f)+0.5
    return E
end

#振り子のシミュレーション
q=Quaternion(1,0,0,0)
ω=[0,0,0]
E=energy(q,ω,(r0,m,g))
ans=[(q,ω,E)]
