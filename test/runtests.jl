using Monolith
using Test
using StaticArrays
using Images
using CSV
using DataFrames
using Glob

import Monolith: Plane, Ray, Rcub, dimensions, create_rcub,
       ray_plane, propagate_ray_rcub,
       points_rcub, vectors_spherical2,
       ggun, imagerun, 
       load_images

th=37.2f0
dd = 50.0f0
X0 = Plane(SVector{3,Float32}(1,0,0), 0.0)
X1 = Plane(SVector{3,Float32}(1,0,0), -dd)
Y0 = Plane(SVector{3,Float32}(0,1,0), 0.0)
Y1 = Plane(SVector{3,Float32}(0,1,0), -dd)
Z0 = Plane(SVector{3,Float32}(0,0,1), 0.0)
Z1 = Plane(SVector{3,Float32}(0,0,1), -th)
mcr = Rcub(X0,X1,Y0,Y1,Z0,Z1)
drc = dimensions(mcr)

rz = Ray(SVector{3,Float64}(10,10,10), SVector{3,Float64}(0,0,1))

#@testset "loadimages" begin
#    load_images("./")
#end

@testset "Rcub" begin
    @test drc.x0 == 0.0
    @test drc.x1 == dd
    @test drc.y0 == 0.0
    @test drc.y1 == dd
    @test drc.z0 == 0.0
    @test drc.z1 == th

    rcbx = create_rcub(drc.x1, drc.y1, drc.z1)
    dimensions(rcbx) == drc
end

@testset "Propagation" begin
    @test ray_plane(rz, mcr.X0) ≈ 1e+7
	@test ray_plane(rz, mcr.X1)≈ 1e+7
	@test ray_plane(rz, mcr.Y0)≈ 1e+7
	@test ray_plane(rz, mcr.Y1)≈ 1e+7

    @test ray_plane(rz, mcr.Z0) == -10.0f0
    @test ray_plane(rz, mcr.Z1) == th - 10.0f0
    @test propagate_ray_rcub(mcr, rz) == th - 10.0f0
end

@testset "generators" begin
    
    drc = dimensions(mcr)
    xpts = points_rcub(mcr, 2)
    size(xpts) == (2,4)
    @test typeof.(xpts.vx)[1] == Float32

    length(vectors_spherical2()) == 3
    typeof(vectors_spherical2()[1])

    xpts =  vectors_spherical2(2) 
    size(xpts) == (2,4)
    @test typeof.(xpts.vx)[1] == Float32
end

@testset "ggun" begin

    lblx, imgx = ggun(drc, ng=500, xsipm=6.0f0, pde=1.0f0, rpos=true)

    @test size(lblx) == (1,4)
    @test lblx.ng[1] == 500
    @test typeof(lblx.xg[1]) == Float32

    @test typeof(imgx) == Matrix{Float32}
    @test size(imgx) == (8,8)
end

@testset "imagerun" begin
    imagerun(1, 2, drc; ipath="./", lpath="./",
             ng=100, xsipm=6.0f0, pde=1.0f0, rpos=true, seed=12344)
    
    rimg1 = Float32.(load("img1.png"))
    rimg2 = Float32.(load("img2.png")) 

    rlbl1 =DataFrame(CSV.File("lbl1.csv"))
    rlbl2 =DataFrame(CSV.File("lbl2.csv"))
    #rlbl1 = Arrow.Table("lbl1")
    #rlbl2 = Arrow.Table("lbl2")

    @test size(rimg1) == (8,8)
    @test size(rimg2) == (8,8)
    @test typeof(rlbl1.ng) == typeof(rlbl2.ng)
    @test typeof(rlbl1.xg) == typeof(rlbl2.xg)


    
end


