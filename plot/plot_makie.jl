using Makie
using FileIO, Colors

function mp_earth()
    earth = try
        load(download("https://svs.gsfc.nasa.gov/vis/a000000/a002900/a002915/bluemarble-2048.png"))
    catch e
        @warn("Downloadinging the earth failed. Using random image, so this test will fail! (error: $e)")
        rand(RGBAf0, 100, 100) # don't error test when e.g. offline
    end
    m = GLNormalUVMesh(Sphere(Point3f0(0), 1f10), 60)
    mesh(m, color = earth, shading = false)
end

function mp_basic_coordinate(r=[0,0,0])
    arrows!([r[1],r[1],r[1]],[r[2],r[2],r[2]],[r[3],r[3],r[3]],[1,0,0],[0,1,0],[0,0,1],normalize=true,linewidth=5,arrowsize=0.1)
end

function mp_show_coodinate(r=[0,0,0],x=[1,0,0],y=[0,1,0],z=[0,0,1])
    
end

function mp_r(r)
    arrows!([0],[0],[0],r[1],r[2],r[3])
end


#Makie.save("/Users/takuya/Desktop/julia_plots/fig.png")


