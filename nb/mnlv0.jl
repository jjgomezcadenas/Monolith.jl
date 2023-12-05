### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ ff4e37ce-d80a-4c79-845f-70eea4f045e2
begin
	using LaTeXStrings
	using StaticArrays
	using LinearAlgebra
	using PlutoUI
	using Printf
	using GLMakie
	using Test
	using Images
	using Random
	using DataFrames
	using DataFramesMeta
end

# ╔═╡ 8c3504ba-8fa0-11ee-175e-5d6c8f47f55b
"""
Defines a ray (e.g, a segment that propagates in straight line)
r = e + t * d, 
where e in the initial point and d is the direction vector
"""
struct Ray
    e::SVector{3,Float64}
    d::SVector{3,Float64}
end

# ╔═╡ a3e9871a-2dfb-426b-b4df-f3af1ef9fd1e
"""
Propagates ray r along the distance defined by t
"""
ray(r::Ray, t::Float64) =r.e + t * r.d

# ╔═╡ 1e7d910c-35fd-4c63-a3fe-2df8c8a540af
"""
Defines a Plane

A x + B y + C  Z + D=0

where norm([A,B,C]) = 1
"""
struct Plane
    P::SVector{3,Float64} # [A, B, C]
    D::Float64
end

# ╔═╡ 601f8388-f233-405e-9d26-6a654c9a40be
md"""
## Example : Intersection of a ray with a plane
- Define a plane [1,0,0, -7] (that is a plane where x = 7)
- Define a ray with origin [2,3,4] and direction [0.577, 0.577, 0.577]
"""

# ╔═╡ b670b8d2-fcb2-4254-9039-81088b6409ef
p7 = Plane(SVector{3,Float64}(1,0,0), -7.0)

# ╔═╡ 6d10f2b2-93ce-40ae-a77f-b366faed7e05
r0 = Ray(SVector{3,Float64}(2,3,4), SVector{3,Float64}(0.577,0.577,0.577))

# ╔═╡ 28bb42f2-e119-4e01-99f9-b729e39caf78
vd = dot(r0.d, p7.P)

# ╔═╡ 57f77743-8583-42ba-8208-ba6d826560bf
md"""
- If vd = 0, the ray is parallel to the plane.
- for vd < 0 Plane points into the ray
- for vd > 0 Plane points away from the ray 
"""

# ╔═╡ 7b3d56c3-ed8b-47c2-8a87-99bbfe2fa695
v0 = -(dot(r0.e, p7.P) + p7.D)

# ╔═╡ 71197aba-f840-4122-ac71-2d2d358a9cb1
md"""
- Now compute t = v0/vd
- if t < 0 the ray goes away from the plane and no intersection occurs
"""

# ╔═╡ 5c9572af-5cdc-4929-83c0-0fc8cb71e03a
t = v0/vd

# ╔═╡ 52dd13b7-5754-473a-b233-218cd9022395
ri = ray(r0,t)

# ╔═╡ 109f13ae-c5a5-48be-8a91-c853ebe8228b
function ray_plane(r::Ray, p::Plane, eps=1e-7, large=1e+7)
	vd = dot(r.d, p.P)
	v0 = -(dot(r.e, p.P) + p.D)
	if abs(vd) < eps
		return large
	else
		return v0/vd
	end
end

# ╔═╡ 3f854fa1-2c70-486d-b374-e40233365063
"""
- Define a rectancular cuboid (Rcub) in terms of 6 planes: (Z0, Zt), (X0, Xd), (Y0, Yd), where:
    - t is the crystal's thickness (e.g, 2 X0~2 x 18.6 for CsI)
    - d is the crystal transverse dimensions (e.g, 50 mm)
- And the planes defining the crystal are:
    - (X0, X1), located at origin and dx: (1,0,0,0) and (1,0,0,-dx)
    - (Y0, Y1), located at origin and dy: (0,1,0,0) and (0,1,0,-dy)
    - (Z0, Z1), located at origin and dz: (0,0,1,0) and (0,0,1,-dz)
"""
struct Rcub
    X0::Plane
    X1::Plane
	Y0::Plane
    Y1::Plane
	Z0::Plane
    Z1::Plane
end

# ╔═╡ 877295d5-c8b3-4bac-9987-f40d9558749b
@bind th NumberField(0.0:0.1:50.0; default=18.6*2)

# ╔═╡ 9d30ccad-7dcf-4c19-84a0-63fc2eb45712
@bind dd NumberField(0.0:0.1:60.0; default=50.0)

# ╔═╡ 8f3a7d2b-364e-4ad5-9dfb-ff0d665283d5
begin
	X0 = Plane(SVector{3,Float64}(1,0,0), 0.0)
	X1 = Plane(SVector{3,Float64}(1,0,0), -dd)
	Y0 = Plane(SVector{3,Float64}(0,1,0), 0.0)
	Y1 = Plane(SVector{3,Float64}(0,1,0), -dd)
	Z0 = Plane(SVector{3,Float64}(0,0,1), 0.0)
	Z1 = Plane(SVector{3,Float64}(0,0,1), -th)
	mcr = Rcub(X0,X1,Y0,Y1,Z0,Z1)
end

# ╔═╡ 6154b124-3d1f-48c2-881d-16560f9db51b
function dimensions(r::Rcub)
	(x0=-X0.D,x1=-X1.D,y0=-Y0.D,y1=-Y1.D,z0=-Z0.D,z1=-Z1.D)
end

# ╔═╡ fa548b66-925e-473c-bafb-91d6a7575c42
dimensions(mcr)

# ╔═╡ 8f2a1377-760a-4749-a94c-f2234325508f
md"""
### Test
- Ray (10,10,10) (0,0,1) along z. Should give "infinite" distance (no crossing) with planes X,Y, t < 0 for Z0 and t >0 for Z1. 
"""

# ╔═╡ 70d3f914-9377-4b94-9b71-68a89f580583
rz = Ray(SVector{3,Float64}(10,10,10), SVector{3,Float64}(0,0,1))

# ╔═╡ dcd72051-2182-41cd-ac76-116493815945
mcr.X0

# ╔═╡ 2b28509d-3701-4a0e-879b-cd2bfe110a31
begin
	@test ray_plane(rz, mcr.X0) ≈ 1e+7
	@test ray_plane(rz, mcr.X1)≈ 1e+7
	@test ray_plane(rz, mcr.Y0)≈ 1e+7
	@test ray_plane(rz, mcr.Y1)≈ 1e+7
end

# ╔═╡ 82dc1c9a-a8d9-4a30-9136-c2729dd755c8


# ╔═╡ 34c4ea3a-7c2e-410c-a1fb-06d2395ac2cf
ray_plane(rz, mcr.Z0)

# ╔═╡ dddf8b8d-2e12-47bc-9d68-666cbbf4435d
ray_plane(rz, mcr.Z1)

# ╔═╡ 8c6f6857-d569-420a-9208-6edf7795b429
md" - And intersection with the plane is: "

# ╔═╡ 834d7280-c9e7-42f8-af10-5412d4c51fb8
ray(rz,ray_plane(rz, mcr.Z1))

# ╔═╡ 6f53d6b7-cc04-4aee-939b-b2231c4c02de
md"""
- Ray (10,10,10) (1,0,0) along X. Should give "infinite" distance (no crossing) with planes Z,Y, t < 0 for X0 and t >0 for X1. 
"""

# ╔═╡ 6875721a-efe3-4a06-9cfc-d2dacc0da7e1
rx = Ray(SVector{3,Float64}(10,10,10), SVector{3,Float64}(1,0,0))

# ╔═╡ dbc07fa9-203e-4e58-8be9-17ed759ea009
begin
	@test ray_plane(rx, mcr.Z0) ≈ 1e+7
	@test ray_plane(rx, mcr.Z1)≈ 1e+7
	@test ray_plane(rx, mcr.Y0)≈ 1e+7
	@test ray_plane(rx, mcr.Y1)≈ 1e+7
end

# ╔═╡ 13fb4601-2daa-412a-bb3a-ec8bce03625d
ray_plane(rx, mcr.X0)

# ╔═╡ 28a4f317-3655-4e7a-ae3e-46e02d271643
ray_plane(rx, mcr.X1)

# ╔═╡ 525067f3-c567-47cf-bacd-92dd218e7c82
ray(rx,ray_plane(rx, mcr.X1))

# ╔═╡ 4ec9f4ed-6016-4018-a80c-99476bc2e2ce
md"""
- Ray (10,10,10) (0,1,0) along Y. Should give "infinite" distance (no crossing) with planes X,Z, t < 0 for Y0 and t >0 for Y1. 
"""

# ╔═╡ 5965b3f0-8242-46c4-b1a6-dc55faf6d0df
ry = Ray(SVector{3,Float64}(10,10,10), SVector{3,Float64}(0,1,0))

# ╔═╡ 4801bcac-d4dc-4e1c-8885-f539412c1de3
begin
	@test ray_plane(ry, mcr.Z0) ≈ 1e+7
	@test ray_plane(ry, mcr.Z1)≈ 1e+7
	@test ray_plane(ry, mcr.X0)≈ 1e+7
	@test ray_plane(ry, mcr.X1)≈ 1e+7
end

# ╔═╡ 136cad33-dcb9-431e-a1b2-b1cb52ca7a5c
ray_plane(ry, mcr.Y0)

# ╔═╡ cb8ce6b1-241d-4c61-acce-bfc04d09ae4c
ray_plane(ry, mcr.Y1)

# ╔═╡ e98392f1-2172-494f-a2a6-7eb30aa59238
ray(ry,ray_plane(ry, mcr.Y1))

# ╔═╡ e6c169b2-8127-44da-99f5-09be0765bce0
md"""
## Drawing with Makie
"""

# ╔═╡ f12d640c-c9b3-49aa-b9a8-d4e3457cbe3f
md"""
## The meshgrid concept

- Consider the figure below. Suppose that we want to define a two-dimensional area formed by the coordinate points x = [0, 1, 2], y = [0, 1, 2]. The coordinates of the resulting 9 points are represented in the figure below. 
"""

# ╔═╡ 95bd0543-08d5-4dc7-88f1-8ab0f2fd0759
mgrid = load("tutorial_numpy_function_meshgrid_02.jpeg")

# ╔═╡ 09905d7e-08aa-4bba-af76-6a3958b1c31f
md"""
- In julia is very simple to obtaind the matrix of tuples representing those point.
- First we define the function f(x,y) = (x,y)
- Then we apply it to the vector represeting x (a column vector of dimensions 3x1) times the adjoint of y (a raw vector of dimension 1x3) to obtain a 3x3 matrix
"""

# ╔═╡ d30e5401-fcda-45ee-9e08-b08e0f7fa9ea
begin 
	xs = [0, 1, 2]
	ys = [0, 1, 2]	
	f(x, y) = (x, y)
	MC = f.(xs,ys')
end

# ╔═╡ 92c00933-4677-450b-ad83-e99000a231f6
md"""
- Notice that in Julia (like in Fortran but unlike C and phython) c the first array index is a row and the second is a column (Julia has column-major arrays). Thus in the matrix of pairs above the first row describes the first column of the figure. 
"""

# ╔═╡ f62adfa7-ed44-4ea8-8755-3d852c7e26fd
md"""
- We now want to obtain two arrays, each of size (3, 3) that include the x and y coordinates (separately) of the 9 points in question. This is also very easy.

ones(3) $\times$ xs' => $\begin{bmatrix}
           1 \\
           1 \\
           1
         \end{bmatrix} (0\,1\,2) \times  = \begin{bmatrix}
           0 \,1 \,2\\
           0 \,1 \,2 \\
           0\, 1 \, 2
         \end{bmatrix}$

 ys $\times$ ones(3)' => $\begin{bmatrix}
           0 \\
           1 \\
           2
         \end{bmatrix} (1\,1\,1) \times  = \begin{bmatrix}
           0 \,0 \,0\\
           1 \,1 \,1 \\
           2\, 2 \, 2
         \end{bmatrix}$
"""

# ╔═╡ df661cdb-cb28-427c-babd-8bfc028fb34b
ones(3) * xs'

# ╔═╡ 3ef1b529-16e9-4ced-96da-d1cc2604e06b
ys * ones(3)'

# ╔═╡ 598416ee-9367-420f-931c-58321d73195a
ones(3)* ones(3)'

# ╔═╡ 4338e485-1f6a-467d-a52a-8184cead0362
md"""
- Notice that the matrices read now by row.

- This is exactly what the numpy.meshgrid function does: it accepts as input the coordinates that define the hyperplane segment (it can be two dimensions or any other number of them) and returns arrays with the coordinates of those points

- Given two vectors, jmgrid returns the julia equivalent to the meshgrid function.

- If transp = true it returns the exact equivalent, which is the matrix transposed of the usual matrix in julia. 
"""

# ╔═╡ 4ca021ad-3fa7-41bc-bb2f-b08a122b7c9f
"""
Given two vectors, jmgrid returns the julia equivalent to the meshgrid function.
If transp = true it returns the exact equivalent, which is the matrix transposed of
the usual matrix in julia. In Julia, contrary to Python, the first array index is a row and the second is a column (Julia has column-major arrays)
"""
function jmgrid(p0::Vector{<:Number}, p1::Vector{<:Number} ; transp=true)
	x = first.(Iterators.product(p0, p1))
	y = last.(Iterators.product(p0, p1))
	if transp
		collect(transpose(x)), collect(transpose(y))
	else
		x,y
	end
end

# ╔═╡ 7bb895a8-7b5d-423b-a4ed-b4caabb24980
jmgrid(xs,ys)

# ╔═╡ 312a78d3-ad4b-4d5f-bd26-61ea17c522e2
jmgrid(xs,ys, transp =false)

# ╔═╡ 3ec2a9e5-4c85-4686-aa51-b4cbb41c6f60
md"""

### Draw the cube defining surfaces
"""

# ╔═╡ c9fa4aaf-0e26-4434-8fb5-90d560d48978
dimensions(mcr)

# ╔═╡ db62d492-99e6-4285-8058-1f308c4842f7
md"""
### Drawing the Rcub
- The idea is to draw the six surfaces that define the cube. A surface can be drawn by using the jmgrid function that returns a meshgrid (two matrices) that define this particular surface. 

- For example the surfaces XY1 and XY2 correspond to planes located at Z = 0 and Z=dz. Those planes are defined with a mesh grid in (XY) which goes between (0,dx,0,dy), then the Z position is a third matrix filled with zeros for the first plane and with the value of dz for the second plane. 
"""

# ╔═╡ 473056d9-407a-4521-9cd6-8c385220fa17
function draw_rcube(mcr::Rcub, pdf, rdf, np=100)
	dm = dimensions(mcr)

	tx = collect(LinRange(dm.x0, dm.x1, np))
	ty = collect(LinRange(dm.y0, dm.y1, np)) 
	tz = collect(LinRange(dm.z0, dm.z1, np)) 
	
	XY1, XY2 = jmgrid(tx, ty)
	XZ1, XZ2 = jmgrid(tx, tz)
	YZ1, YZ2 = jmgrid(ty, tz)
	
	O = ones(np)* ones(np)'
	#GLMakie.surface(X, Y, Z)
	SX0 = (XY1, XY2, O .* dm.z0)
	SX1 = (XY1, XY2, O .* dm.z1)

	SY0 = (XZ1, O .* dm.x0, XZ2)
	SY1 = (XZ1, O .* dm.x1, XZ2)

	SZ0 = (O .* dm.y0, YZ1, YZ2)
	SZ1 = (O .* dm.y1, YZ1, YZ2)

	#SX0, SX1, SY0, SY1, SZ0, SZ1
	scene = Scene()
	cam3d!(scene)
	sfx0 = GLMakie.surface!(scene, SX0..., transparency = true, fxaa = true, alpha = 0.2, overdraw=false, colormap=:reds)
	sfx1 = GLMakie.surface!(scene, SX1..., transparency = true, fxaa = true, alpha = 0.2, overdraw=false, colormap=:reds)
	sfy0 = GLMakie.surface!(scene, SY0..., transparency = true, fxaa = true, alpha = 0.2, overdraw=false, colormap=:blues)
	sfy1 = GLMakie.surface!(scene, SY1..., transparency = true, fxaa = true, alpha = 0.2, overdraw=false, colormap=:blues)
	sfz0 = GLMakie.surface!(scene, SZ0..., transparency = true, fxaa = true, alpha = 0.2, overdraw=false, colormap=:greens)
	sfz1 = GLMakie.surface!(scene, SZ1..., transparency = true, fxaa = true, alpha = 0.2, overdraw=false, colormap=:greens)

	
	scatter!(scene, pdf.vx, pdf.vy, pdf.vz)
	scatter!(scene, rdf.vx, rdf.vy, rdf.vz)

	for (i, r) in enumerate(eachrow(pdf)) 
		lx = [pdf[i,"vx"], rdf[i,"vx"]]
		ly = [pdf[i,"vy"], rdf[i,"vy"]]
		lz = [pdf[i,"vz"], rdf[i,"vz"]]
		
		lines!(scene, lx, ly, lz)
	end
	center!(scene)
	scene
end

# ╔═╡ 4a7cea23-68d5-4c44-9ffa-14d7dc174975
md"""
### Generating gammas inside crystal.

- (X,Y,Z) generated randomly between crystal limits

"""

# ╔═╡ 2ca5aaaf-5a42-4873-84be-aadd48df9133
dm = dimensions(mcr)

# ╔═╡ 94d881e2-86d8-430e-a93c-1daa8da64ef5
Pg = [rand(dm.x0:0.1:dm.x1), rand(dm.y0:0.1:dm.y1), rand(dm.z0:0.1:dm.z1)]

# ╔═╡ 679f59c2-b4af-4031-9886-7f07278796bf
"""
Generate random points inside the cuboid


"""
function points_rcub(mcr::Rcub, npoints::Integer, seed=12345)
    Random.seed!(seed)
	dm = dimensions(mcr)

	vx = rand(dm.x0:0.1:dm.x1, npoints)
	vy = rand(dm.y0:0.1:dm.y1, npoints)
	vz = rand(dm.z0:0.1:dm.z1,npoints)
	DataFrame(ID=1:npoints, vx=vx, vy=vy, vz=vz)
end

# ╔═╡ da90427d-076b-439c-9b82-e5c131c67c4d
pdf = points_rcub(mcr, 10)

# ╔═╡ d25f814f-f80a-4fb4-b1ca-eb2c049810ed
pdf.vx

# ╔═╡ b6529a7a-de00-4926-ac74-34727bdb2622


# ╔═╡ 64dbc2a9-27f1-45a7-a873-de54ebfb7990
"""
Generate three standard normally distributed numbers 
and normalize the vector.

We have to be careful in the case that the vector 
has a norm close to zero, in which we must worry about 
floating point precision by dividing by a very small number. 
This is the reason for the while loop.

"""
function vectors_spherical2(npoints::Integer, eps=1e-5, seed=12345)
    Random.seed!(seed)
    
	vx = zeros(npoints)
	vy = zeros(npoints)
	vz = zeros(npoints)

	for i = 1:npoints
		v = vectors_spherical2(eps, seed, false)
		vx[i] = v[1]
		vy[i] = v[2]
		vz[i] = v[3]
	end
	 DataFrame(ID=1:npoints, vx=vx, vy=vy, vz=vz)
end


# ╔═╡ 62b671aa-fd63-4f4f-a292-34ec32fb96c5
function vectors_spherical2(eps=1e-5, seed=12345, actseed=true)
    if actseed 
		Random.seed!(seed)
	end
  	
	v = zeros(3)
	while norm(v) < eps
		x = randn()  # random standard normal
		y = randn()
		z = randn()
		v = [x, y, z]
	end
    
    v / norm(v)
		
end

# ╔═╡ 95ada1a7-8b61-4370-b92c-46187d23dcb8
norm(vectors_spherical2())

# ╔═╡ e4087b65-b235-4c39-a9d5-a8260720b0c3
vdf = vectors_spherical2(10)

# ╔═╡ 4b8c29cb-5973-45a2-aca0-483d30d056d6
function propagate_ray_rcub(mcr, rz, large=1e+7)
	dst = zeros(6)
	dst[1] = ray_plane(rz, mcr.X0)
	dst[2] = ray_plane(rz, mcr.X1)
	dst[3] = ray_plane(rz, mcr.Y0)
	dst[4] = ray_plane(rz, mcr.Y1)
	dst[5] = ray_plane(rz, mcr.Z0)
    dst[6] = ray_plane(rz, mcr.Z1)

	I = []
	for i in 1:6
		#@info "i" i
		#@info "dst[i]" dst[i]
		if dst[i] > 0 && dst[i] < large
			append!(I,i)
		end
	end
	#@info "I" I
	xmin = large
	for ii in I
		if dst[ii] < xmin
			xmin = dst[ii]
		end
	end
	xmin
end

# ╔═╡ 04eef2f5-adb0-4810-904e-fbf16781ed57
propagate_ray_rcub(mcr, rz)

# ╔═╡ ba436b72-8aaf-4ace-accd-e02d7340dec2
for (i, r) in enumerate(eachrow(pdf)) 
    	
		@info "pdf" r.vx, r.vy, r.vz
		@info "vdf" vdf[i,"vx"], vdf[i,"vy"],vdf[i,"vz"]
		ry = Ray(SVector{3,Float64}(r.vx,r.vy,r.vz), 
			     SVector{3,Float64}(vdf[i,"vx"], vdf[i,"vy"],vdf[i,"vz"]))
		@info "ray" ry
		t = propagate_ray_rcub(mcr,ry)
		@info "t" t
		P = ray(ry,t)
		@info "P =", P 
	end

# ╔═╡ f9af1088-ae40-428f-bee6-e8815281d43b
md"""

## Putting it all together

- Draw the crystal.
- Generate and draw gammas
- Transport gammas
"""

# ╔═╡ fa56ef37-12a9-4cac-9e77-7fa42f4354b9
function create_rcub(dx::Float64,dy::Float64,dz::Float64)
	X0 = Plane(SVector{3,Float64}(1,0,0), -0.0)
	X1 = Plane(SVector{3,Float64}(1,0,0), -dx)
	Y0 = Plane(SVector{3,Float64}(0,1,0), -0.0)
	Y1 = Plane(SVector{3,Float64}(0,1,0), -dy)
	Z0 = Plane(SVector{3,Float64}(0,0,1), -0.0)
	Z1 = Plane(SVector{3,Float64}(0,0,1), -dz)
	Rcub(X0,X1,Y0,Y1,Z0,Z1)
end

# ╔═╡ e3675491-680c-4f66-a110-34a537d6ccb7
function sicrys(dx::Float64,dy::Float64,dz::Float64, ng=10; plot=false)
	mcr = create_rcub(dx,dy,dz) # create cube
	pdf = points_rcub(mcr, ng) # positions
	vdf =vectors_spherical2(ng) # direction vectors

	dm = dimensions(mcr)
	@info dm
	#propagate gammas


	px = zeros(ng)
	py = zeros(ng)
	pz = zeros(ng)

	nz = 0
	for (i, r) in enumerate(eachrow(pdf)) 
    	
		#@info "pdf" r.vx, r.vy, r.vz
		#@info "vdf" vdf[i,"vx"], vdf[i,"vy"],vdf[i,"vz"]
		ry = Ray(SVector{3,Float64}(r.vx,r.vy,r.vz), 
			     SVector{3,Float64}(vdf[i,"vx"], vdf[i,"vy"],vdf[i,"vz"]))
		#@info "ray" ry
		t = propagate_ray_rcub(mcr,ry)
		#@info "t" t
		P = ray(ry,t)
		
		#@info "P =", P 
		px[i] = P[1]
		py[i] = P[2]
		pz[i] = P[3]

		if px[i] > dm.x1
			px[i] = dm.x1
		end
		if py[i] > dm.y1
			py[i] = dm.y1
		end
		if pz[i] > dm.z1
			pz[i] = dm.z1
		end

		if px[i] < dm.x0
			px[i] = dm.x0
		end
		if py[i] < dm.y0
			py[i] = dm.y0
		end
		if pz[i] < dm.z0
			pz[i] = dm.z0
		end

	end

	
	rdf = DataFrame(ID=1:ng, vx=px, vy=py, vz=pz)
	
	if plot
		draw_rcube(mcr, pdf, rdf)
	else
		pdf, rdf
	end
end

# ╔═╡ 71589e75-c6b8-4619-9e7c-6266cc0de894


# ╔═╡ 41280d6f-77f4-45ec-9b8e-04811d175873
ipdf, idf = sicrys(50.0,50.0,18.6*2,10000, plot=false)

# ╔═╡ 469427f9-4d28-4f65-b0f8-b944131513c7
18.6 * 2

# ╔═╡ 16d70580-5f2b-48c3-bdf0-dc04f152d9d4
dfpc = @chain idf begin
		@rsubset :vz > 36.9 
		end

# ╔═╡ 32f806b4-6774-4ea5-acbb-7fc9dbdde795
hist(idf.vz; bins = 50, normalization = :none)

# ╔═╡ 0bf147dc-ef24-456c-b2d4-a2d332efef16
scatter(dfpc.vx, dfpc.vy)



# ╔═╡ e488957a-5be7-4264-a0c9-5bc0a0595393
md"""
## Shooting one point at the time, shooting ng gammas and recording the SiPM signal
"""

# ╔═╡ 83f57912-971c-4543-9298-255d92b30009
function crgun(dx::Float64,dy::Float64,dz::Float64; ng=10^3, xsipm=3)
	function getxyz(dm::NamedTuple)
		rand(dm.x0:0.1:dm.x1), rand(dm.y0:0.1:dm.y1), rand(dm.z0:0.1:dm.z1)
	end
	
	function fill_point!(P::SVector{3,Float64}, px::Vector{Float64}, 
	                     py::Vector{Float64},   pz::Vector{Float64}, 
		                  i::Integer,           dm::NamedTuple)
		px[i] = P[1]
		py[i] = P[2]
		pz[i] = P[3]

		if px[i] > dm.x1
			px[i] = dm.x1
		end
		if py[i] > dm.y1
			py[i] = dm.y1
		end
		if pz[i] > dm.z1
			pz[i] = dm.z1
		end

		if px[i] < dm.x0
			px[i] = dm.x0
		end
		if py[i] < dm.y0
			py[i] = dm.y0
		end
		if pz[i] < dm.z0
			pz[i] = dm.z0
		end

	end
	#mcr = create_rcub(dx,dy,dz) # create cube
	#dm = dimensions(mcr)

	@info "Crystal with dimensions" dm

	x,y,z = getxyz(dm) # generate position

	@info "Shooting from position" x y z

	@info "Shooting number of gammas ng =" ng
	
	vdf =vectors_spherical2(ng) # direction vectors
	
	#propagate gammas

	px = zeros(ng)
	py = zeros(ng)
	pz = zeros(ng)

	nz = 0
	for (i, _) in enumerate(eachrow(vdf)) 
    	
		ry = Ray(SVector{3,Float64}(x,y,z), 
			     SVector{3,Float64}(vdf[i,"vx"], vdf[i,"vy"],vdf[i,"vz"]))	
		t = propagate_ray_rcub(mcr,ry)
		fill_point!(ray(ry,t), px, py, pz, i, dm)		
	end

	rdf = DataFrame(ID=1:ng, vx=px, vy=py, vz=pz)
	dfpc = @chain rdf begin
		@rsubset :vz == dz 
		@select!  $(Not(:vz))
	end
	
	dfpc
end

# ╔═╡ bf6c1ca1-57b0-4c5d-93e8-fbc687eb5804
ixdf = crgun(50.0,50.0,18.6*2,ng=10^3, xsipm=3)

# ╔═╡ 9ad1c431-5d59-483f-b2b1-979b6ee1f0a8
scatter(ixdf.vx, ixdf.vy)

# ╔═╡ e017e966-77ba-4feb-bfea-4c8a8e291288
function crgun2(dm::NamedTuple; x::Float64, y::Float64, z::Float64, 
                ng=10^3, xsipm=3,seed=12345, eps=1e-5)
	
	function fill_point!(P::SVector{3,Float64}, px::Vector{Float64}, 
	                     py::Vector{Float64},   pz::Vector{Float64}, 
		                  i::Integer,           dm::NamedTuple)
		px[i] = P[1]
		py[i] = P[2]
		pz[i] = P[3]

		if px[i] > dm.x1
			px[i] = dm.x1
		end
		if py[i] > dm.y1
			py[i] = dm.y1
		end
		if pz[i] > dm.z1
			pz[i] = dm.z1
		end

		if px[i] < dm.x0
			px[i] = dm.x0
		end
		if py[i] < dm.y0
			py[i] = dm.y0
		end
		if pz[i] < dm.z0
			pz[i] = dm.z0
		end

	end
	#mcr = create_rcub(dx,dy,dz) # create cube
	#dm = dimensions(mcr)

	
	@info "Crystal with dimensions" dm

	@info "Shooting from position" x y z

	@info "Shooting number of gammas ng =" ng
	
	vdf =vectors_spherical2(ng, eps, seed) # direction vectors
	
	#propagate gammas

	px = zeros(ng)
	py = zeros(ng)
	pz = zeros(ng)

	nz = 0
	for (i, _) in enumerate(eachrow(vdf)) 
    	
		ry = Ray(SVector{3,Float64}(x,y,z), 
			     SVector{3,Float64}(vdf[i,"vx"], vdf[i,"vy"],vdf[i,"vz"]))	
		t = propagate_ray_rcub(mcr,ry)
		fill_point!(ray(ry,t), px, py, pz, i, dm)		
	end

	rdf = DataFrame(ID=1:ng, vx=px, vy=py, vz=pz)
	dfpc = @chain rdf begin
		@rsubset :vz == dm.z1 
		@select!  $(Not(:vz))
	end
	
	dfpc
end

# ╔═╡ afd5c875-738b-420b-b55c-4a9b3addf174
@bind yg NumberField(0.0:0.1:50.0; default=25.0)

# ╔═╡ 684a1ddd-2e10-4433-bb46-06831f780f14
@bind zg NumberField(0.0:0.1:37.2; default=18.0)

# ╔═╡ 820f125a-1373-4c96-ab47-c84e40339997
@bind gMeV NumberField(0:10^5; default=500000)

# ╔═╡ e4be4572-4442-4239-817a-97ee8980b4bd
ixdf2 = crgun2(dm, x=xg, y=yg, z=zg, ng=gMeV, xsipm=3, seed=1236);

# ╔═╡ 36f4d5d1-720b-4c3f-80e1-9dcaf9119153
size(ixdf2)

# ╔═╡ fb2169de-3f60-4903-a111-3b82b9b9511f
scatter(ixdf2.vx, ixdf2.vy)

# ╔═╡ 284c2bbc-1c95-4091-b153-d391f64691a7
@bind xsipm NumberField(0:0.1:9; default=3.0)

# ╔═╡ 3fec99f3-749c-4671-aed7-68b919f1e00f
begin
	nbx = Int(floor(dm.x1 /xsipm))
	nby = Int(floor(dm.y1 /xsipm))
	nbox = nbx * nby
	md"""
	- number of boxes along x = $nbx
	- number of boxes along y = $nby
	- total number of boxes = $nbox
	"""
end

# ╔═╡ 84e08531-4211-4bd1-a306-91b20cfa77c6
function indxy(xsipm::Float64)
	function indx(x::Vector{Float64}) 
		xx = Int.(floor.(x  ./xsipm ) )
		replace(xx, 0 =>1) 
	end
	indx
end

# ╔═╡ 2653bcb4-12e2-4740-acc6-2eadaab97634
xindx = indxy(3.0)

# ╔═╡ ed2e0da8-5563-4d7b-9484-66cc3780d921
xindx([0.0, 50.0])

# ╔═╡ 52660d04-2c4d-4e04-9d53-f89ced64572b
dfs = @chain ixdf2 begin 
    transform(_, [:vx]  => xindx => :ix)
	transform(_, [:vy]  => xindx => :iy)
	end

# ╔═╡ 2e46a12f-67fc-41de-8614-8313587b94f5
# 5. Group gammas in boxes of the size of the SiPm
	dfbx = unique(combine(groupby(dfs, [:ix, :iy])) do xypdf
				nb = nrow(xypdf) 
				xbox = sum(xypdf.vx) / nb
				ybox = sum(xypdf.vy) / nb
				DataFrame(
					      ix =xypdf.ix,
					      iy =xypdf.iy,
				          ebox = nb, xbox=xbox, ybox=ybox)
    end)
		 


# ╔═╡ 996ba008-ebbf-4c2b-b6ae-424b21f89141
dfbx.ix

# ╔═╡ 9019f016-c456-46e3-931d-3e66a27d3dab
dfbx.iy

# ╔═╡ b8c16539-f667-456b-9e34-79748d362425
iar = argmax(dfbx.ebox)

# ╔═╡ 212923cd-6620-42b5-998f-5c4f8288336d
dfbx.ebox

# ╔═╡ 177a3b9c-a755-48a2-bd8e-46377c3caddb
dfbx[iar,:]

# ╔═╡ 9af594cb-5510-4b50-82c9-0f8b18346028
exy = zeros(nbx,nby)

# ╔═╡ e75445ef-b03d-4095-9fd8-921a5d4c2f50
xr = zeros(nbx)

# ╔═╡ 6b84ed45-1ece-4236-88c2-5a77d7d4646c
yr = zeros(nby)

# ╔═╡ 55ce56e2-dd00-4708-b8e2-4b1d98da7025
for r in eachrow(dfbx) 
	exy[r.ix,r.iy] = r.ebox
	xr[r.ix] = r.xbox
	yr[r.iy] = r.ybox
end

# ╔═╡ 05f85a49-9aa7-42ab-ac12-acf1caafeb41
imgxy =exy ./ maximum(dfbx.ebox)

# ╔═╡ ebed9cad-c03e-40af-a426-d3d3d2e90b5b
length(xr)

# ╔═╡ 71b55a56-53a3-4e5b-b93b-3fa6450020a4
length(yr)

# ╔═╡ e4176907-04ad-40c4-b71a-b6fec6633f14
size(exy)

# ╔═╡ d90df5bc-4b4a-492f-8e91-343946223b62
surface(xr, yr, exy)

# ╔═╡ 30e49368-bbb7-4cae-a507-3c679a07fda1
image(imgxy)

# ╔═╡ cf7b07a7-b53c-400f-b88e-c4c50aad2300
@bind xg NumberField(0.0:0.1:50.0; default=25.0)

# ╔═╡ 2c47ccff-4d63-429b-9f4b-7683e53d0796
# ╠═╡ disabled = true
#=╠═╡
xg = rand(0:0.1:50.0,100)
  ╠═╡ =#

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
GLMakie = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[compat]
DataFrames = "~1.6.1"
DataFramesMeta = "~0.14.1"
GLMakie = "~0.9.1"
Images = "~0.25.3"
LaTeXStrings = "~1.3.1"
PlutoUI = "~0.7.54"
StaticArrays = "~1.7.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.4"
manifest_format = "2.0"
project_hash = "97b7d324d7f188b63a5042529387cbf8b4003cd5"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractLattices]]
git-tree-sha1 = "222ee9e50b98f51b5d78feb93dd928880df35f06"
uuid = "398f06c4-4d28-53ec-89ca-5b2656b7603d"
version = "0.3.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "793501dcd3fa7ce8d375a2c878dca2296232686e"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "02f731463748db57cc2ebfbd9fbc9ce8280d3433"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.7.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "247efbccf92448be332d154d6ca56b9fcdd93c31"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.6.1"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "0da671c730d79b8f9a88a391556ec695ea921040"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.0.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "32abd86e3c2025db5172aa182b982debed519834"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.1"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "e0af648f0692ec1691b5d094b8724ba1346281cf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.18.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "05f9816a77231b07e634ab8715ba50e5249d6f76"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.5"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "f9d7112bfff8a19a3a4ea4e03a8e6a91fe8456bf"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "Reexport"]
git-tree-sha1 = "6970958074cd09727b9200685b8631b034c0eb16"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.14.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["DataStructures", "EnumX", "ExactPredicates", "Random", "SimpleGraphs"]
git-tree-sha1 = "26eb8e2331b55735c3d305d949aabd7363f07ba7"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "0.8.11"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "a6c00f894f24460379cb7136633cef54ac9f6f4a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.103"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ErrorfreeArithmetic]]
git-tree-sha1 = "d6863c556f1142a061532e79f611aa46be201686"
uuid = "90fa49ef-747e-5e6f-a989-263ba693cf1a"
version = "0.5.2"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArraysCore"]
git-tree-sha1 = "499b1ca78f6180c8f8bdf1cabde2d39120229e5c"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.6"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.Extents]]
git-tree-sha1 = "2140cd04483da90b2da7f99b2add0750504fc39c"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FastRounding]]
deps = ["ErrorfreeArithmetic", "LinearAlgebra"]
git-tree-sha1 = "6344aa18f654196be82e62816935225b3b9abe44"
uuid = "fa42c844-2597-5d31-933b-ebd51ab2693f"
version = "0.3.1"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "299dc33549f68299137e51e6d49a13b5b1da9673"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.1"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "28e4e9c4b7b162398ec8004bdabe9a90c78c122d"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.8.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "c6e4a1fbe73b31a3dea94b1da449503b8830c306"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.21.1"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "50351f83f95282cf903e968d7c6e8d44a5f83d0b"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW]]
deps = ["GLFW_jll"]
git-tree-sha1 = "35dbc482f0967d8dceaa7ce007d16f9064072166"
uuid = "f7f18e0c-5ee9-5ccd-a5bf-e8befd85ed98"
version = "3.4.1"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GLMakie]]
deps = ["ColorTypes", "Colors", "FileIO", "FixedPointNumbers", "FreeTypeAbstraction", "GLFW", "GeometryBasics", "LinearAlgebra", "Makie", "Markdown", "MeshIO", "ModernGL", "Observables", "PrecompileTools", "Printf", "ShaderAbstractions", "StaticArrays"]
git-tree-sha1 = "875190f74d7694020a10c6c8f9377d2d4f8e7a5c"
uuid = "e9467ef8-e4e7-5192-8a1a-b1aee30e663a"
version = "0.9.1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "d53480c0793b13341c40199190f92c611aa2e93c"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.2"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "424a5a6ce7c5d97cca7bcc4eac551b97294c54af"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.9"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "899050ace26649433ef1af25bc17a815b3db52b7"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.9.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "af13a277efd8a6e716d79ef635d5342ccb75be61"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.10.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageContrastAdjustment]]
deps = ["ImageBase", "ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "eb3d4365a10e3f3ecb3b115e9d12db131d28a386"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.12"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "08b0e6354b21ef5dd5e49026028e41831401aca8"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.17"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "PrecompileTools", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "3447781d4c80dbe6d71d239f7cfb1f8049d4c84f"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.6"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "bca20b2f5d00c4fbc192c3212da8fa79f4688009"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.7"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils"]
git-tree-sha1 = "b0b765ff0b4c3ee20ce6740d843be8dfce48487c"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.3.0"

[[deps.ImageMagick_jll]]
deps = ["JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1c0a2295cca535fabaf2029062912591e9b61987"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.10-12+3"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.ImageMorphology]]
deps = ["ImageCore", "LinearAlgebra", "Requires", "TiledIteration"]
git-tree-sha1 = "e7c68ab3df4a75511ba33fc5d8d9098007b579a8"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.3.2"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "LazyModules", "OffsetArrays", "PrecompileTools", "Statistics"]
git-tree-sha1 = "783b70725ed326340adf225be4889906c96b8fd1"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.7"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "44664eea5408828c03e5addb84fa4f916132fc26"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.8.1"

[[deps.ImageShow]]
deps = ["Base64", "ColorSchemes", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "3b5344bcdbdc11ad58f3b1956709b5b9345355de"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.8"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "8717482f4a2108c9358e5c3ca903d3a6113badc9"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.9.5"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "5fa9f92e1e2918d9d1243b1131abe623cdf98be7"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.25.3"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "ea8031dea4aff6bd41f1df8f2fdfb25b33626381"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.4"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "be8e690c3973443bec584db3346ddc904d4884eb"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.5"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ad37c091f7d7daf900963171600d7c1c5c3ede32"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "EnumX", "FastRounding", "LinearAlgebra", "Markdown", "Random", "RecipesBase", "RoundingEmulator", "SetRounding", "StaticArrays"]
git-tree-sha1 = "f59e639916283c1d2e106d2b00910b50f4dab76c"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.21.2"

[[deps.IntervalSets]]
deps = ["Dates", "Random"]
git-tree-sha1 = "3d8866c029dd6b16e69e0d4a939c4dfcb98fac47"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.8"
weakdeps = ["Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "PrecompileTools", "Printf", "Reexport", "Requires", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "9bbb5130d3b4fa52846546bca4791ecbdfb52730"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.38"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "d65930fa2bc96b07d7691c652d701dcbe7d9cf0b"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "90442c50e202a5cdf21a7899c66b240fdef14035"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.7"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "3a994404d3f6709610701c7dabfc03fed87a81f8"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.1"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "Permutations", "Primes", "SimplePolynomials"]
git-tree-sha1 = "8889b8aa6821c0ee73828a2139314b4d596e7dbc"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.2.2"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "eb006abbd7041c28e0d16260e50a24f8f9104913"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Setfield", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "7f9d0e484fb9090d383283d8528e3ffab987f661"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.20.1"

[[deps.MakieCore]]
deps = ["Observables", "REPL"]
git-tree-sha1 = "e81e6f1e8a0e96bf2bf267e4bf7f94608bf09b5c"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.7.1"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "96ca8a313eb6437db5ffe946c457a401bbb8ce1d"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MeshIO]]
deps = ["ColorTypes", "FileIO", "GeometryBasics", "Printf"]
git-tree-sha1 = "8be09d84a2d597c7c0c34d7d604c039c9763e48c"
uuid = "7269a6da-0436-5bbc-96c2-40638cbb6118"
version = "0.4.10"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "1130dbe1d5276cb656f6e1094ce97466ed700e5a"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.ModernGL]]
deps = ["Libdl"]
git-tree-sha1 = "b76ea40b5c0f45790ae09492712dd326208c28b2"
uuid = "66fc600b-dfda-50eb-8b99-91cfa97b1301"
version = "1.1.7"

[[deps.Mods]]
git-tree-sha1 = "61be59e4daffff43a8cec04b5e0dc773cbb5db3a"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "1.3.3"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.Multisets]]
git-tree-sha1 = "8d852646862c96e226367ad10c8af56099b4047e"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "2ac17d29c523ce1cd38e27785a7d23024853a4bb"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.10"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cc6e1927ac521b659af340e0ca45828a3ffc748f"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.12+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "01f85d9269b13fedc61e63cc72ee2213565f7a72"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.8"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "4e5be6bb265d33669f98eb55d2a57addd1eeb72c"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.30"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "eed372b0fa15624273a9cdb188b1b88476e6a233"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.2"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "a935806434c9d4c506ba941871b327b96d41f2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.0"

[[deps.Permutations]]
deps = ["Combinatorics", "LinearAlgebra", "Random"]
git-tree-sha1 = "c7745750b8a829bc6039b7f1f0981bcda526a946"
uuid = "2ae35dd2-176d-5d53-8349-f30d82d94d4f"
version = "0.4.19"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "bd7c69c7f7173097e7b5e1be07cee2b8b7447f51"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.54"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "a9c7a523d5ed375be3983db190f6a5874ae9286d"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.6"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "1d05623b5952aed1307bf8b43bec8b8d1ef94b6e"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.5"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9ebcd48c498668c7fa0e97a9cae873fbee7bfee1"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.1"

[[deps.Quaternions]]
deps = ["LinearAlgebra", "Random", "RealDot"]
git-tree-sha1 = "9a46862d248ea548e340e30e2894118749dc7f51"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.7.5"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RealDot]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "9f0a1b71baaf7650f4fa8a1d168c7fb6ee41f0c9"
uuid = "c1ae055f-0cd5-4b69-90a6-9a35b1a98df9"
version = "0.1.0"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.RingLists]]
deps = ["Random"]
git-tree-sha1 = "f39da63aa6d2d88e0c1bd20ed6a3ff9ea7171ada"
uuid = "286e9d63-9694-5540-9e3c-4e6708fa07b2"
version = "0.2.8"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays"]
git-tree-sha1 = "792d8fd4ad770b6d517a13ebb8dadfcac79405b8"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.6.1"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "0e7508ff27ba32f26cd459474ca2ede1bc10991f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SetRounding]]
git-tree-sha1 = "d7a25e439d07a17b7cdf97eecee504c50fedf5f6"
uuid = "3cc68bcd-71a2-5612-b932-767ffbe40ab0"
version = "0.2.1"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "db0219befe4507878b1a90e07820fed3e62c289d"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.0"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleGraphs]]
deps = ["AbstractLattices", "Combinatorics", "DataStructures", "IterTools", "LightXML", "LinearAlgebra", "LinearAlgebraX", "Optim", "Primes", "Random", "RingLists", "SimplePartitions", "SimplePolynomials", "SimpleRandom", "SparseArrays", "Statistics"]
git-tree-sha1 = "f65caa24a622f985cc341de81d3f9744435d0d0f"
uuid = "55797a34-41de-5266-9ec1-32ac4eb504d3"
version = "0.8.6"

[[deps.SimplePartitions]]
deps = ["AbstractLattices", "DataStructures", "Permutations"]
git-tree-sha1 = "e9330391d04241eafdc358713b48396619c83bcb"
uuid = "ec83eff0-a5b5-5643-ae32-5cbf6eedec9d"
version = "0.3.1"

[[deps.SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "e14e1a7063179a90c2981faf3c8cbf12aacccdaf"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.16"

[[deps.SimpleRandom]]
deps = ["Distributions", "LinearAlgebra", "Random"]
git-tree-sha1 = "3a6fb395e37afab81aeea85bae48a4db5cd7244a"
uuid = "a6525b86-64cd-54fa-8f65-62fc48bdc0e8"
version = "0.3.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays"]
git-tree-sha1 = "4b33e0e081a825dbfaf314decf58fa47e53d6acb"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.4.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "5165dfb9fd131cf0c6957a3a7605dede376e7b63"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.0"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableHashTraits]]
deps = ["Compat", "SHA", "Tables", "TupleTools"]
git-tree-sha1 = "d29023a76780bb8a3f2273b29153fd00828cb73f"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "1.1.1"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "5ef59aea6f18c25168842bded46b16662141ab87"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.7.0"
weakdeps = ["Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["Adapt", "ConstructionBase", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "0a3db38e4cce3c54fe7a71f831cd7b6194a54213"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.16"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "34cc045dd0aaa59b8bbe86c644679bc57f1d5bd0"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.8"

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

[[deps.TranscodingStreams]]
git-tree-sha1 = "1fbeaaca45801b4ba17c251dd8603ef24801dd84"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.2"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "155515ed4c4236db30049ac1495e2969cc06be9d"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.4.3"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "5f24e158cf4cee437052371455fe361f526da062"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.6"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "da69178aacc095066bad1f69d2f59a60a1dd8ad1"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.0+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╠═ff4e37ce-d80a-4c79-845f-70eea4f045e2
# ╠═8c3504ba-8fa0-11ee-175e-5d6c8f47f55b
# ╠═a3e9871a-2dfb-426b-b4df-f3af1ef9fd1e
# ╠═1e7d910c-35fd-4c63-a3fe-2df8c8a540af
# ╠═601f8388-f233-405e-9d26-6a654c9a40be
# ╠═b670b8d2-fcb2-4254-9039-81088b6409ef
# ╠═6d10f2b2-93ce-40ae-a77f-b366faed7e05
# ╠═28bb42f2-e119-4e01-99f9-b729e39caf78
# ╠═57f77743-8583-42ba-8208-ba6d826560bf
# ╠═7b3d56c3-ed8b-47c2-8a87-99bbfe2fa695
# ╠═71197aba-f840-4122-ac71-2d2d358a9cb1
# ╠═5c9572af-5cdc-4929-83c0-0fc8cb71e03a
# ╠═52dd13b7-5754-473a-b233-218cd9022395
# ╠═109f13ae-c5a5-48be-8a91-c853ebe8228b
# ╠═3f854fa1-2c70-486d-b374-e40233365063
# ╠═877295d5-c8b3-4bac-9987-f40d9558749b
# ╠═9d30ccad-7dcf-4c19-84a0-63fc2eb45712
# ╠═8f3a7d2b-364e-4ad5-9dfb-ff0d665283d5
# ╠═6154b124-3d1f-48c2-881d-16560f9db51b
# ╠═fa548b66-925e-473c-bafb-91d6a7575c42
# ╠═8f2a1377-760a-4749-a94c-f2234325508f
# ╠═70d3f914-9377-4b94-9b71-68a89f580583
# ╠═dcd72051-2182-41cd-ac76-116493815945
# ╠═2b28509d-3701-4a0e-879b-cd2bfe110a31
# ╠═82dc1c9a-a8d9-4a30-9136-c2729dd755c8
# ╠═04eef2f5-adb0-4810-904e-fbf16781ed57
# ╠═34c4ea3a-7c2e-410c-a1fb-06d2395ac2cf
# ╠═dddf8b8d-2e12-47bc-9d68-666cbbf4435d
# ╠═8c6f6857-d569-420a-9208-6edf7795b429
# ╠═834d7280-c9e7-42f8-af10-5412d4c51fb8
# ╠═6f53d6b7-cc04-4aee-939b-b2231c4c02de
# ╠═dbc07fa9-203e-4e58-8be9-17ed759ea009
# ╠═13fb4601-2daa-412a-bb3a-ec8bce03625d
# ╠═28a4f317-3655-4e7a-ae3e-46e02d271643
# ╠═525067f3-c567-47cf-bacd-92dd218e7c82
# ╠═6875721a-efe3-4a06-9cfc-d2dacc0da7e1
# ╠═4ec9f4ed-6016-4018-a80c-99476bc2e2ce
# ╠═5965b3f0-8242-46c4-b1a6-dc55faf6d0df
# ╠═4801bcac-d4dc-4e1c-8885-f539412c1de3
# ╠═136cad33-dcb9-431e-a1b2-b1cb52ca7a5c
# ╠═cb8ce6b1-241d-4c61-acce-bfc04d09ae4c
# ╠═e98392f1-2172-494f-a2a6-7eb30aa59238
# ╠═e6c169b2-8127-44da-99f5-09be0765bce0
# ╠═f12d640c-c9b3-49aa-b9a8-d4e3457cbe3f
# ╠═95bd0543-08d5-4dc7-88f1-8ab0f2fd0759
# ╠═09905d7e-08aa-4bba-af76-6a3958b1c31f
# ╠═d30e5401-fcda-45ee-9e08-b08e0f7fa9ea
# ╠═92c00933-4677-450b-ad83-e99000a231f6
# ╠═f62adfa7-ed44-4ea8-8755-3d852c7e26fd
# ╠═df661cdb-cb28-427c-babd-8bfc028fb34b
# ╠═3ef1b529-16e9-4ced-96da-d1cc2604e06b
# ╠═598416ee-9367-420f-931c-58321d73195a
# ╠═4338e485-1f6a-467d-a52a-8184cead0362
# ╠═4ca021ad-3fa7-41bc-bb2f-b08a122b7c9f
# ╠═7bb895a8-7b5d-423b-a4ed-b4caabb24980
# ╠═312a78d3-ad4b-4d5f-bd26-61ea17c522e2
# ╠═3ec2a9e5-4c85-4686-aa51-b4cbb41c6f60
# ╠═c9fa4aaf-0e26-4434-8fb5-90d560d48978
# ╠═db62d492-99e6-4285-8058-1f308c4842f7
# ╠═473056d9-407a-4521-9cd6-8c385220fa17
# ╠═4a7cea23-68d5-4c44-9ffa-14d7dc174975
# ╠═2c47ccff-4d63-429b-9f4b-7683e53d0796
# ╠═2ca5aaaf-5a42-4873-84be-aadd48df9133
# ╠═94d881e2-86d8-430e-a93c-1daa8da64ef5
# ╠═679f59c2-b4af-4031-9886-7f07278796bf
# ╠═da90427d-076b-439c-9b82-e5c131c67c4d
# ╠═d25f814f-f80a-4fb4-b1ca-eb2c049810ed
# ╠═b6529a7a-de00-4926-ac74-34727bdb2622
# ╠═64dbc2a9-27f1-45a7-a873-de54ebfb7990
# ╠═62b671aa-fd63-4f4f-a292-34ec32fb96c5
# ╠═95ada1a7-8b61-4370-b92c-46187d23dcb8
# ╠═e4087b65-b235-4c39-a9d5-a8260720b0c3
# ╠═4b8c29cb-5973-45a2-aca0-483d30d056d6
# ╠═ba436b72-8aaf-4ace-accd-e02d7340dec2
# ╠═f9af1088-ae40-428f-bee6-e8815281d43b
# ╠═fa56ef37-12a9-4cac-9e77-7fa42f4354b9
# ╠═e3675491-680c-4f66-a110-34a537d6ccb7
# ╠═71589e75-c6b8-4619-9e7c-6266cc0de894
# ╠═41280d6f-77f4-45ec-9b8e-04811d175873
# ╠═469427f9-4d28-4f65-b0f8-b944131513c7
# ╠═16d70580-5f2b-48c3-bdf0-dc04f152d9d4
# ╠═32f806b4-6774-4ea5-acbb-7fc9dbdde795
# ╠═0bf147dc-ef24-456c-b2d4-a2d332efef16
# ╠═e488957a-5be7-4264-a0c9-5bc0a0595393
# ╠═83f57912-971c-4543-9298-255d92b30009
# ╠═bf6c1ca1-57b0-4c5d-93e8-fbc687eb5804
# ╠═9ad1c431-5d59-483f-b2b1-979b6ee1f0a8
# ╠═e017e966-77ba-4feb-bfea-4c8a8e291288
# ╠═cf7b07a7-b53c-400f-b88e-c4c50aad2300
# ╠═afd5c875-738b-420b-b55c-4a9b3addf174
# ╠═684a1ddd-2e10-4433-bb46-06831f780f14
# ╠═820f125a-1373-4c96-ab47-c84e40339997
# ╠═e4be4572-4442-4239-817a-97ee8980b4bd
# ╠═36f4d5d1-720b-4c3f-80e1-9dcaf9119153
# ╠═fb2169de-3f60-4903-a111-3b82b9b9511f
# ╠═284c2bbc-1c95-4091-b153-d391f64691a7
# ╠═3fec99f3-749c-4671-aed7-68b919f1e00f
# ╠═84e08531-4211-4bd1-a306-91b20cfa77c6
# ╠═2653bcb4-12e2-4740-acc6-2eadaab97634
# ╠═ed2e0da8-5563-4d7b-9484-66cc3780d921
# ╠═52660d04-2c4d-4e04-9d53-f89ced64572b
# ╠═2e46a12f-67fc-41de-8614-8313587b94f5
# ╠═996ba008-ebbf-4c2b-b6ae-424b21f89141
# ╠═9019f016-c456-46e3-931d-3e66a27d3dab
# ╠═b8c16539-f667-456b-9e34-79748d362425
# ╠═212923cd-6620-42b5-998f-5c4f8288336d
# ╠═177a3b9c-a755-48a2-bd8e-46377c3caddb
# ╠═9af594cb-5510-4b50-82c9-0f8b18346028
# ╠═e75445ef-b03d-4095-9fd8-921a5d4c2f50
# ╠═6b84ed45-1ece-4236-88c2-5a77d7d4646c
# ╠═55ce56e2-dd00-4708-b8e2-4b1d98da7025
# ╠═05f85a49-9aa7-42ab-ac12-acf1caafeb41
# ╠═ebed9cad-c03e-40af-a426-d3d3d2e90b5b
# ╠═71b55a56-53a3-4e5b-b93b-3fa6450020a4
# ╠═e4176907-04ad-40c4-b71a-b6fec6633f14
# ╠═d90df5bc-4b4a-492f-8e91-343946223b62
# ╠═30e49368-bbb7-4cae-a507-3c679a07fda1
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
