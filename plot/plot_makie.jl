using Makie
using FileIO, Colors

"""
mp_init()

makie plotの初期化

"""
function mp_init()
    scene = Scene()
end


"""
mp_earth(earth_radius=1)

地球を表示する関数

# Arguments
- （earth_radius:地球の半径）
"""
function mp_earth(;earth_radius=1)
    earth = try
        load(download("https://svs.gsfc.nasa.gov/vis/a000000/a002900/a002915/bluemarble-2048.png"))
    catch e
        @warn("Downloadinging the earth failed. Using random image, so this test will fail! (error: $e)")
        rand(RGBAf0, 100, 100) # don't error test when e.g. offline
    end
    m = GLNormalUVMesh(Sphere(Point3f0(0), earth_radius), 60)
    mesh(m, color = earth, shading = false)
end


"""
mp_basic_coordinate(scene,r=[0,0,0])

基底を示す矢印を表示・重ねて表示

# Arguments
- myscene：上書き対象となるscene
- （r:座標系の原点）
- （arrow_size:矢印の長さ）

"""
function mp_basic_coordinate(myscene; r=[0,0,0],arrow_size=1)
    arrows!(myscene,[r[1],r[1],r[1]],[r[2],r[2],r[2]],[r[3],r[3],r[3]],[arrow_size,0,0],[0,arrow_size,0],[0,0,arrow_size],normalize=false,linewidth=5,arrowsize=arrow_size/10)
end


function mp_r(myscene, r)
    arrows!(myscene, [0],[0],[0],r[1],r[2],r[3])
end

function mp_sat(myscene,r)
    mp_basic_coordinate(myscene)
    arrows!(myscene,[0],[0],[0],[r[1]],[r[2]],[r[3]],normalize=true,linewidth=5,arrowsize=0.1)
end


"""
mp_orbit(r)

軌道の表示
# Arguments
- myscene：上書き対象となるscene
- r_log：軌道の記録

"""
function mp_orbit(myscene, r_log)
    mp_basic_coordinate(myscene,arrow_size=10000)
    for i = 1:size(r_log)[1]
        arrows!(myscene,[0],[0],[0],[r_log[i,1]],[r_log[i,2]],[r_log[i,3]],normalize=false,linewidth=5,arrowsize=0.1)
    end

    Base.prompt("終了します")
end


"""
mp_orbit_normal(r)

軌道の表示
# Arguments
- myscene：上書き対象となるscene
- r_log：軌道の記録

"""
function mp_orbit_normal(myscene, r_log)
    mp_basic_coordinate(myscene,arrow_size=1)
    for i = 1:size(r_log)[1]
        arrows!(myscene,[0],[0],[0],[r_log[i,1]],[r_log[i,2]],[r_log[i,3]],normalize=true,linewidth=1,arrowsize=0.01)
    end

    Base.prompt("終了します")
end


"""
mp_sat_env(r)

軌道面座標系上での各パラメータの表示
# Arguments
- myscene：上書き対象となるscene
- vecs：表示したいベクトル（正規化される）

"""
function mp_sat_env(myscene, vecs,vec1...)
    mp_basic_coordinate(myscene,arrow_size=1)
    for i = 1:size(vecs)[1]
        arrows!(myscene,[0],[0],[0],[vecs[i,1]],[vecs[i,2]],[vecs[i,3]],normalize=true,linewidth=1,arrowsize=0.05)
    end

    Base.prompt("終了します")
end


function mp_delete(myscene,nums)
    for i = 1:size(nums)
        delete!(myscene, myscene[i])
    end
end

#Makie.save("/Users/takuya/Desktop/julia_plots/fig.png")


