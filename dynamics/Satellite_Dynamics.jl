#衛星ダイナミクス
include("DynamicEnvironmentField.jl")
using DifferentialEquations
using SatelliteToolbox
using Plots;gr()
using Quaternions
using DEF
using LinearAlgebra
#パラメータ
ω0=[0;0;0]
q0=quat(1,1,1,0)
u0=[q0,ω0]
tspan=(0.0,100,0)#時間
Ix=1.0
Iy=1.0
Iz=1.0
I=
[Ix 0 0
;0 Iy 0
;0 0 Iz]
invI=inv(I)
#常微分方程式の定義
du(u,invI,t)=[(1/2)u[1]*u[2],invI*DEF.Torque(u)]
plob=ODEProblem(du,u0,tspan)
#状微分方程式の求解
sol=solve(plob)
