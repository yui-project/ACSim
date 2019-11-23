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
x0=[q0,ω0]
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
dx(x,t)=[(1/2)x[1]*x[2],invI*DEF.Torque(x)]
plob=ODEProblem(dx,x0,tspan)
#状微分方程式の求解
sol=solve(plob)

