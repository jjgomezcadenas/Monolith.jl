using Monolith
using Test
using StaticArrays
import Monolith: Plane, Rcub, dimensions, points_rcub

function test_rcub()

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
    
    @test drc.x0 == 0.0
    @test drc.x1 == dd
    @test drc.y0 == 0.0
    @test drc.y1 == dd
    @test drc.z0 == 0.0
    @test drc.z1 == th
end

function test_points_rcub()
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
    xpts = points_rcub(mcr, 2)
    size(xpts) == (2,4)
end
