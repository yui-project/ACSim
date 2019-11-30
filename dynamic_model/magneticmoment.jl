using LinearAlgebra  

function f(I1,I2,I3,ml,B)
    #Ii:電流の大きさを持った、コイル電流右ねじ向きのベクトル
    #ml:残留磁気モーメント
    #B:地磁気
 t1=cross(I1,B)
 t2=cross(I2,B)
 t3=cross(I3,B) 
 Tl=cross(ml,B)

 (0.00000081)(22/7)*(t1+t2+t3)+Tl
end

println(
    f(
        [1,2,3],[4,5,6],[5,6,2],[4,3,2],[4,5,6]
        )
        )

