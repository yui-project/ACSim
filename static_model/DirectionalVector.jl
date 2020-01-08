using LinearAlgebra
function f(r,p)
    # r: position of Satellite 
    # p: target

    d =  p - r
    # directional vector
    m = norm(d) 
    #norm of d

    d_hat=(1/m)d
    # directional vector (normalized)
   
   end