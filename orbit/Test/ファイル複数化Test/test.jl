include("Func.jl")
using Plots
gr()


a = 1.0
b = 2.0
c = [1. ,2. ,3.]
Data = [1., 2., 3, 4., 5., 6.]
Num =0.
u = 0.0:0.1:1.0
v = 8. .*u

function mul(x,y)
    z = x*y
    return z
end


##別ファイルで定義した関数が使えるか
println(plus(a, b))
### →使える

##"."で要素ごとの演算ができる性質は維持されるのか
println(plus.(c, b))
### →維持される

##同じファイルで複数の関数を定義した際にそれぞれを使えるか
println(minus(a,b))
### →使える

##元ファイルで定義した変数を関数ファイルで使えるか
Print()
### →使える

##関数内で、元ファイルの変数を上書きできるか
println(Num)
OverWrite()
println(Num)
### →できない

##元ファイルで定義した関数を関数ファイルで使えるか
println(Pow(b, 5))
### →使える

##元ファイルで導入したモジュールが関数ファイルで使えるか
Pl(u,v)
### →使えない(多分?)


