
#using Plots
#gr()
module MyModule
    
    export plus
    function plus(x,y)
        z = x + y
        return z
    end

    function minus(x,y)
        z = x - y
        return z
    end

    function Print()
        for i=1:6
            println(Data[i])
        end
    end

    function OverWrite()
        Num = 1.
        println(Num)
    end

    function Pow(x,y)
        z = 1
        for i=1:y
            z = mul(z,x)
        end
        return z
    end

    function Pl(u,v)
        println(u)
        println(v)
        plot(u,v)
        savefig("test3.png")
        tles = read_tle("ISS_TLE.txt") #TLE読み込み
        println(tles.ω)
    end
end
