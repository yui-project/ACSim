include("differencial.jl")

function RK4(F,dt)
    f=dif(F)
    F=F[1]
    k1=dt*f
    k2=dt*(F+k1)
    k3=dt*(F+0.5(k1+k2))
    k4=dt*(F+k1+k3)
    return F+(k1+2k2+2k3+k4)/6
end
