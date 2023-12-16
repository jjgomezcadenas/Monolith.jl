### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ a298b03b-c7c4-4114-9984-d37e7df8d4c8
using Pkg; Pkg.activate("/Users/jjgomezcadenas/Projects/JRT")

# ╔═╡ 5158389a-fca1-11ed-05a1-c7345c9bd354
begin
	#using Revise
	using PlutoUI
	using CSV
	using DataFrames
	using DataFramesMeta
	#using Plots
	#using Printf
	using Markdown
	using InteractiveUtils
	using Statistics
	using Random
	using Distributions
	using PDMats
	using Chain
	using LinearAlgebra
	using Images
	using GLMakie
	#using Makie
	import GeometryBasics
	using Colors
	using Roots
	using Unitful
	using UnitfulEquivalences
end

# ╔═╡ e24e5710-04b9-49c3-b3f6-2b136a1744eb
import Unitful:
    nm, μm, mm, cm, m, km,
    mg, g, kg,
    ps, ns, μs, ms, s, minute, hr, d, yr, Hz, kHz, MHz, GHz,
    eV,
    μJ, mJ, J,
	μW, mW, W,
	A, N, mol, mmol, V, L, M,
    Bq, mBq, μBq 

# ╔═╡ e36751ff-a521-4529-95c8-6bbfc3314e66
begin
	A_BI214_316Ti   =    1.0    * mBq/kg
	A_TL208_316Ti   =    0.4   * mBq/kg
	A_BI214_CU_LIM  =   12      * μBq/kg
	A_TL208_CU_LIM  =    1.4    * μBq/kg
end

# ╔═╡ df3dde7a-9437-409f-bd51-db975f74b78c
PlutoUI.TableOfContents(title="JRT draft", indent=true)


# ╔═╡ fc869138-1cf7-463b-aca5-729741f8b1ff
GLMakie.activate!()

# ╔═╡ f230deec-f4c9-429a-89f5-92d8c08e8472
md"""
### Defining a cylinder using GeometryBasics
"""

# ╔═╡ 7371f833-a7a8-482c-925b-730d8378c22e
gcyl = GeometryBasics.Cylinder(GeometryBasics.Point{3, Float64}(1,2,3), GeometryBasics.Point{3, Float64}(2,3,4), 1.0)

# ╔═╡ 2c803aef-9a52-4108-a252-e63a95929c67
cyL =GeometryBasics.mesh(GeometryBasics.Tesselation(gcyl, 100))

# ╔═╡ c7583c88-4cfc-49e1-b61e-c299801be450
md"""
## Drawing a cylinder
- A cylinder of radius r, and height h, having z-axis as symmetry axis, is a stack of circles of radius r:

$\begin{align}
x= r \times \cos(u) \\
y=r \times sin(u) \\
z=v 
\end{align}$

where:
$u \in [0,2\pi]$
$v \in [0,h]$

- This allows to plot the cylinder as any parameterized surface.
- For example:
	- Define u as a range beteen 0 and 2π
	- Define v as a range between 0 and 10.0 (the length of the cylinder)
"""

# ╔═╡ c48c5f53-2632-4cfc-b07b-7fde9f1d8eb6
md"""
## Testing functions
"""

# ╔═╡ de0fb8fd-c757-4769-a2da-d874add0a553
md"""
### Check function **in_endcaps**
"""

# ╔═╡ c95dcf93-40ff-4db0-a323-29991f5135dd
md"""
### Demonstrate propagation functions
"""

# ╔═╡ c7200fa5-31eb-458c-80ca-70f968edd9b0
md"""
#### Ray moving along z axis, with no x,y components. It must intersect endcups (either zmax or zmin)
"""

# ╔═╡ a89c0c98-ee2c-4103-8cd3-9be03a2b75ae
md"""
- Function **cylinder\_intersection\_roots** returns the smaller positive root corresponding to the intersection between the ray and the cylinder. Since in this case the ray does not intersect the cylinder (the endcups are not the cylinder surface) the root must be zero 
"""

# ╔═╡ aa445ed2-dc2d-4802-ab12-18505487ceef
md"""
- Function **ray\_intersection\_with_cylinder** will return the recomputed ray
"""

# ╔═╡ 6c2b1465-02ae-4768-a4f6-e8f461339b97
md"""
- Since this is a ray going in the zp direction, t=0, and the recomputed e point in the ray is in zmax and belongs to the upper end-cap. 
"""

# ╔═╡ 26d72ea7-31e5-4c55-af55-7666f7c47f1c
md"""
- We can now check that the point **inendcaps** is indeed in the end-caps
"""

# ╔═╡ b9c450a6-446e-4a5c-93ee-28a45cf96356
md"""
- A ray pointing into negative z should end up in the bottom end-cap
"""

# ╔═╡ 59138131-46c1-4ac8-be0b-42f5ca2e30b8
md"""
#### Ray intersecting the cylinder surface
"""

# ╔═╡ 02a7fd08-440c-4e0c-93c8-d42792642d11
md"""
- The ray is moving with large director cosines so it will intersect the cylinder surface.
"""

# ╔═╡ bca3de2e-74a6-423e-a660-3845c3e9657b
md"""
- We can visualize where the roots will be
"""

# ╔═╡ c97578f2-1d39-497a-9e79-31e041e77b7b
md"""
#### Ray intersecting the cylinder surface much further than zmax
"""

# ╔═╡ 68e9ab79-21a4-48c1-82c3-45d8862e70cf
md"""
### Demonstrate generation of photons
"""

# ╔═╡ b4110420-ef1c-4e1b-868a-68bd9d943fea
md"""
- Points must be normalized
"""

# ╔═╡ 5e866260-a798-4c5f-b407-38d091d0ed07
md"""
- Points must be in a unit sphere
"""

# ╔═╡ be8dd99e-60f5-4d75-959c-cffbd5bfa4b8
md"""
### Demonstrate generation around one point and propagation to cylinder
"""

# ╔═╡ cf2cce55-84fa-408a-b7ff-6748754eea34
begin
	rcx = 70.0/2.0 #cm
	lcx = 100.0 #cm
	zcrx = 5.0 
	pcrx = atan(zcrx, rcx)
end

# ╔═╡ 8a45b62d-cbfb-4f8c-9c61-fad0e3ca53da
md"""
- Generate photons around one point and transport photons to cylinder surfaces
"""

# ╔═╡ 128d5ed6-3a09-427e-b50d-159d0288d1dd
np=100000

# ╔═╡ 8baa74f5-557e-49cb-a6c2-a18be5e3c1d1
# Transforms a vector of values between 0.0 and dx (dy) in 
	# a vector of indices
	function indxy(xsipm::Float64)
		function indx(x::Vector{Float64}) 
			xx = Int.(floor.(x  ./xsipm ) )
			replace(xx, 0 =>1) 
		end
		indx
	end

# ╔═╡ 63bedd95-3e1b-4786-8a14-91db66e503fd
begin
	zindx = indxy(zcrx)
	pindx = indxy(pcrx)
end

# ╔═╡ 31628f89-3dcc-461a-81b2-582420805f91
pcrx

# ╔═╡ 571fe2c8-4b2b-4436-b84b-f59bf9d32b13
(2π/pcrx) * zcrx

# ╔═╡ b6a5d9b3-60d2-44b9-ab35-cb298466b0d4
2π*rcx

# ╔═╡ 35f50d8e-e56a-4e9a-8288-e15b133f065e
#dfs2 = @select(dfs, :gex, :gey, :gez)

# ╔═╡ 5cca07ea-8740-4899-bca7-629bcdb8d46b
begin
ocmax = 2e-3
latt = 9.7e-2 #cm
d_torso = 40 # cm
lgamma = d_torso/2.0
att = exp(-latt*lgamma)
end

# ╔═╡ 5db4cd9d-71bf-4614-b6b6-922f3dcbdbf1
pos = Poisson(1.4 * 0.14)

# ╔═╡ 953f8162-fb20-4062-b185-3271340a135f
pdf(pos, 2)

# ╔═╡ ce9da661-5f2f-4e5d-825e-5d029c274066
md"""
# Types
"""

# ╔═╡ 24b1b8ef-10ee-4d96-8121-827234233873
"""
Simple representation of a physical material

# Fields
- `ρ::typeof(1.0g/cm^3)`          : density (ρ)
- `μovrρ::typeof(1.0cm^2/g)`      : μ/ρ
- `μ::Unitful.Length^-1`          : attenuation coefficient (μ) 
- `λ::Unitful.Length`             : attenuation length 

"""
struct PhysicalMaterial
    ρ::typeof(1.0g/cm^3) 
    μovrρ::typeof(1.0cm^2/g) 
	μ::typeof(1.0/cm)
	λ::Unitful.Length

    function PhysicalMaterial(ρ::typeof(1.0g/cm^3),  μovrρ::typeof(1.0cm^2/g))
        μ = μovrρ*ρ
        λ = 1.0/μ
		new(ρ, μovrρ, μ, λ)
	end
end

# ╔═╡ 6458b9a6-e0b1-40cf-87bd-3e3f01990525
"""
Simple representation of a radioactive material

# Fields
- `m::PhysicalMaterial`          : a radioactive material is a physical material
- `a_bi214::typeof(1.0*Bq/kg)`   : activity of Bi-214
- `a_tl208::typeof(1.0*Bq/kg)`   : activity of Tl-208

"""
struct RadioactiveMaterial
    m::PhysicalMaterial
    a_bi214::typeof(1.0*Bq/kg)
    a_tl208::typeof(1.0*Bq/kg)
    
end

# ╔═╡ 478e01bc-d764-4af7-a6f4-661cdba4820d
function transmittance_at_qbb(rmt::RadioactiveMaterial, d::Float64) 
	l = uconvert(mm, rmt.m.λ)/mm
	exp(-d/l)
end

# ╔═╡ efe04df4-cb56-4916-96f7-296dc37492ac
struct Cylinder
    r   :: Float64
    zmin:: Float64
    zmax:: Float64
    p0::Vector{Float64}
    p1::Vector{Float64}

    function Cylinder(r:: Float64, zmin::Float64, zmax::Float64)
        p0 = [0.0, 0.0, zmin] # point in one endcup
        p1 = [0.0, 0.0, zmax] # point in the other
		new(r, zmin, zmax, p0, p1)
        #mag, v, n1, n2 = unit_vectors_()
        #P, P2, P3 = surfaces_(r, mag, v, n1, n2)
	end
end

# ╔═╡ 3c398d08-ef94-4202-9b87-62aa4fba07ed
cyl =Cylinder(5.0, 0., 10.0) # draws only the barrel 

# ╔═╡ 50f7c64a-4fa6-4705-b2e9-896ef545029f
cyl

# ╔═╡ c7d4775a-f585-4852-a172-337655cb21ed
stdsk = Cylinder(rcx, 0., lcx)

# ╔═╡ f99d7b04-2fc6-48c5-a846-601e2e270691
gp = [0.0,0.0,(stdsk.zmax + stdsk.zmin)/2]

# ╔═╡ 2f01a131-d184-4ea8-97c1-89309716383f
struct PhysicalCylinder
	cyl::Cylinder
	units::Unitful.Length
end

# ╔═╡ 053952a1-0558-447c-a3f5-ae30c2d2f7be
"""
Defines a ray (e.g, a segment that propagates in straight line)
r = e + t * d, 
where e in the initial point and d is the direction vector
"""
struct Ray
    e::Vector{Float64}
    d::Vector{Float64}
	u::Vector{Float64}
	function Ray(e::Vector{Float64}, d::Vector{Float64})
		u = d .- e
		new(e,d,u)
	end
end


# ╔═╡ 9e193e81-0147-46d2-943f-1526dfb26e50
rzp = Ray([1.0, 0.0, 0.0], [0.0, 0.0, 1.0]) # starts in (1,0,0), moves towards z+

# ╔═╡ 5de835fa-284f-4dd0-ba1d-033d4c50da5a
rzn = Ray([1.0, 0.0, 0.0], [0.0, 0.0, -1.0]) # starts in (1,0,0), moves towards z+

# ╔═╡ 40cff386-97de-4e6c-8d22-216302831893
ry = Ray([1.0, 0.0, 0.0], [0.5, 0.5, 1.0])

# ╔═╡ 4d17bf89-d01b-4f1b-9dc7-9dcd00ee1d8a
rz = Ray([1.0, 0.0, 0.0], [0., 0.1, 1.0])

# ╔═╡ 991a026a-bfe2-4bc8-a200-31dd25a881d1
md"""
# Functions
"""

# ╔═╡ 2facda59-dc5a-4951-8d6f-1d121690aad2
"""
Propagates ray r along the distance defined by t
"""
ray(r::Ray, t::Float64) =r.e + t * r.d


# ╔═╡ 539faa01-299b-4c5b-8113-fc21501dd84d
"""
Takes a Cylinder type and returns a GeometryBasics cylinder and mesh. 
"""
function get_cylinder(cl::Cylinder, np=100)
	p0 = GeometryBasics.Point{3, Float64}(cl.p0[1],cl.p0[2],cl.p0[3]) 
	p1 = GeometryBasics.Point{3, Float64}(cl.p1[1],cl.p1[2],cl.p1[3]) 
	gcyl = GeometryBasics.Cylinder(p0, p1, cl.r)
	cyL =GeometryBasics.mesh(GeometryBasics.Tesselation(gcyl, np))
	gcyl, cyL
end

# ╔═╡ ef077643-f5b1-4244-af7e-8111ede962a0
gcyl2, cylm = get_cylinder(cyl) # from cyl get GeometryBasics objects

# ╔═╡ 269e4b09-3505-40fe-8815-2bc847d02a99
"""
Takes a cylinder and returns a Makie drawing, either frame or wireframe
"""
function draw_cylinder(cl::Cylinder, np=100;  resolution = (1000,1000), 
                       col = :blue, clevel=0.5, linewidth = 2.5, wf=true)
	gcyl2, cylm = get_cylinder(cyl,np)

	with_theme(theme_light()) do
		fig = Figure(resolution = (1200,800))
		axs = Axis3(fig[1,1]; aspect=:data, perspectiveness=0.5) 
		if wf == false
    		GLMakie.mesh!(axs, cylm, color = (col, clevel), 
			  transparency = true, shading = false)
		else
    		GLMakie.wireframe!(axs, cylm; color = col, linewidth = linewidth)
		end
		fig
	end
end

# ╔═╡ 81d62300-1a88-4c9c-b050-8acd20857867
fig = draw_cylinder(cyl, wf=true) # draw cylinder as wireframe

# ╔═╡ 5e1c520f-fb11-4698-a981-5f3863c4788c
draw_cylinder(cyl, wf=false) # draw cylinder as volume 

# ╔═╡ c26c8511-afae-4447-9458-d1d21f565911
"""
Equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
"""
cylinder_equation(c::Cylinder, P::Vector{Float64}) = P[1]^2 + P[2]^2 - c.r^2

# ╔═╡ deefebf9-07a3-4cce-b2d5-4b42644b5db7
"""
Length of cylinder
"""
clength(c::Cylinder) = c.zmax - c.zmin



# ╔═╡ 6f6cd50f-d4fd-48d8-8b53-7edf9112c45b
"""
Length of physical cylinder
"""
clength(c::PhysicalCylinder) = clength(c.cyl) * c.units

# ╔═╡ e9d59365-ea79-4e97-b908-cbd528dc0804
"""
Perimeter of cylinder
"""
perimeter(c::Cylinder) = 2π * c.r



# ╔═╡ c66ca1b8-1ae8-4b53-b828-88cac5d8eae9
"""
Perimeter of physical cylinder
"""
perimeter(c::PhysicalCylinder) = perimeter(c.cyl) * c.units

# ╔═╡ da14195e-414e-40ab-aca1-83d2d68b9172
"""
Cylinder area (barrel)
"""
area_barrel(c::Cylinder) = perimeter(c) * clength(c)


# ╔═╡ 0f55b2c3-db7b-4aad-bd8d-c678a5db0b02
"""
Physical Cylinder area (barrel)
"""
area_barrel(c::PhysicalCylinder) = perimeter(c) * clength(c)

# ╔═╡ f41fbfdf-4cee-4b82-b087-e224c70b7f6b
"""
Cylinder area (endcap)
"""
area_endcap(c::Cylinder) = π * c.r^2
        


# ╔═╡ 8da08d4a-f0cf-4b36-b634-b9f8716a3184
"""
Physical Cylinder area (endcap)
"""
area_endcap(c::PhysicalCylinder) = area_endcap(c.cyl) * c.units^2

# ╔═╡ 7f8c81be-0f2d-4214-b27d-a7470534d867
"""
Cylinder area (total)
"""
area(c::Cylinder) = area_barrel(c) + area_endcap(c)
       


# ╔═╡ e1f26c48-0c0c-480c-9574-96b29b65ccf7
"""
Physical Cylinder area (total)
"""
area(c::PhysicalCylinder) = area_barrel(c) + area_endcap(c)

# ╔═╡ 372b5f7a-2919-4c73-8ce4-3ebf8ebb21b9
"""
Cylinder volume 
"""
volume(c::Cylinder) = area_endcap(c) * clength(c)


# ╔═╡ 834027c9-0a39-450f-9996-a1ec56f080a9
"""
Physical Cylinder volume 
"""
volume(c::PhysicalCylinder) = area_endcap(c) * clength(c)

# ╔═╡ 9022e8da-fe82-4631-8e83-549dad785f9a
"""
Mass of a cylinder filled with radioactive material mat
"""
mass(c::PhysicalCylinder, mat::RadioactiveMaterial) = uconvert(g, volume(c) * mat.m.ρ)

# ╔═╡ 8569dc46-642f-4709-9a78-0292f5aa1613
"""
Activity (Bi214) of a cylinder filled with radioactive material mat
"""
a_bi214(c::PhysicalCylinder, mat::RadioactiveMaterial) = uconvert(Bq, mat.a_bi214 * mass(c, mat))


# ╔═╡ 1061c7ff-fbf4-4dd4-bfc3-fb2c38022d87
"""
Activity (Tl208) of a cylinder filled with radioactive material mat
"""
a_tl208(c::PhysicalCylinder, mat::RadioactiveMaterial) = uconvert(Bq, mat.a_tl208 * mass(c, mat))


# ╔═╡ 0abeba90-ba75-4231-9dd8-d35400964c60
 """Returns True if point is in end-caps of cyinder
 This means that the third coordinate of p (p[3], is close to zmin or to zmax)
 """
in_endcaps(c::Cylinder, p::Vector{Float64}) = isapprox(p[3], c.zmin, atol=1e-7) || isapprox(p[3], c.zmax, atol=1e-7) 

# ╔═╡ d06eaf38-f8bb-4732-9b8f-92f2a1053176
in_endcaps(cyl, [0.0, 0.0, 5.0]) # not in endcup

# ╔═╡ 2940f9d9-f2c3-4772-8208-2375370b7e61
in_endcaps(cyl, [0.0, 0.0, 0.0] )  # in bottom endcup

# ╔═╡ a7378f23-9b56-4a5f-8c4b-b51545e69633
in_endcaps(cyl, [0.0, 0.0, 10.0] )  # in top endcup

# ╔═╡ 464a3cf9-c030-4270-bae7-bf4339ed3fc8
 """Returns True if point is in barrel of cyinder
 This means that the point verifies the equation of cylinder
 """
in_barrel(c::Cylinder, p::Vector{Float64}) = isapprox(cylinder_equation(c, p), 0.0, atol= 1.0e-7)

# ╔═╡ 58491933-49ec-4b70-ab1f-2301657ed9e8
"""
Returns true if the point is inside the cylinder
"""
function inside_cylinder(c::Cylinder, p::Vector{Float64})
	s1 =sqrt(p[1]^2 + p[2]^2) <= c.r
	s2 = c.zmin <= p[3] <= c.zmax
	s1 && s2
end

# ╔═╡ 4d53eaca-0cad-4f34-8205-216e9c0a16d9
md"""
To compute the intersection roots we write the ray as:

$r = \mathbf{e} + t \mathbf{d}$

and substitute coordinate by coordinate in the equation of the cylinder:

$F(x,y,z) = x^2 + y^2 - r^2 = 0$

thus:

$(e_1 + t d_1)^2 + (e_2 + t d_2)^2 -r^2 =0$

Then re-arrange

$(d_1^2 + d_2^2)t^2 + 2 t (e_1 d_1 + e_2 d_2) + e_1^2 + e_2^2 - r^2 =0$

So this is a polynomial of the form

$a t^2 + b t + c =0$

with:

$a = (d_1^2 + d_2^2)$

$b = 2 (e_1 d_1 + e_2 d_2)$

$c = e_1^2 + e_2^2 - r^2$

"""

# ╔═╡ 989cec01-548a-40e8-9a82-149dbfa7d367
"""
Computes  the polynomial defined by the intersection roots between a ray and a cylinder. Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
"""
function cylinder_intersection_poly(r::Ray, cl::Cylinder)
    function gf()
		a = r.d[1]^2 + r.d[2]^2
		b = 2 * (r.e[1] * r.d[1] + r.e[2] * r.d[2])
		c = r.e[1]^2 + r.e[2]^2 - cl.r^2

		f(x) =  a * x^2 + b * x + c
	end
	gf()
end

# ╔═╡ 7f4ab967-b5c2-46d9-90c4-ebbbdab0dbfe
ff = cylinder_intersection_poly(ry, cyl)

# ╔═╡ b78b7410-886f-47db-b8a3-b8a7feeec341
lines(-10:0.01:10, ff)

# ╔═╡ c2d3e54d-3ad5-46b6-9d28-be8acf1086b2
find_zeros(ff, -10.0, 10.0)  #finds two roots

# ╔═╡ c6f63bc1-6d4d-49b3-988f-ccaaaae562e1
"""
Computes intersection roots between a ray and a cylinder
Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
"""
function cylinder_intersection_roots(r::Ray, cl::Cylinder, eps::Float64=1e-9)
	
	f = cylinder_intersection_poly(r, cl)
	
	#println("a = ", -2*clength(cyl) * cyl.r, " b = ", 2*clength(cyl) * cyl.r )
	
	rts = find_zeros(f, -100*clength(cyl) * cyl.r, 100*clength(cyl) * cyl.r)
	#println("roots: rts = ", rts)
	if length(rts) == 0
		return 0
	else
		prts = [x for x in rts if x >eps]
		if length(prts) == 0
			return 0
		else
			return minimum(prts)
		end
	end
end

# ╔═╡ da50ac8a-8381-490b-a000-b6ba8b392888
cylinder_intersection_roots(rzp, cyl)

# ╔═╡ 71203991-d5c4-439f-8dd1-fe1b4cf3b875
cylinder_intersection_roots(rzn, cyl)

# ╔═╡ 96af216d-c0e3-4825-a11f-3aaa04c42d0b
cylinder_intersection_roots(ry, cyl) # we keep only the positive root

# ╔═╡ d137c26c-bf75-4899-ac0a-eb8e70c792bb
cylinder_intersection_roots(rz, cyl) 

# ╔═╡ c40ad59c-7687-4bba-b923-fa498295821c
"""
Intersection between a ray and the end-cups of a cylinder
"""
function ray_intersection_with_cylinder_end_caps(r::Ray, cl::Cylinder, t::Float64)
 
    p = ray(r, t)
    if p[3] > cl.zmax
        t = (cl.zmax - r.e[3])/r.d[3]
    else
        t = (cl.zmin - r.e[3])/r.d[3]
	end

    t, ray(r,t)
end

# ╔═╡ 8547b5f7-6b8a-4a9b-b0bf-5b1d8b08c072
 """
 Intersection between a ray and a cylinder
 """
function ray_intersection_with_cylinder(r::Ray, cl::Cylinder)
   
    t = cylinder_intersection_roots(r, cl)

	if t == 0 && r.d[3] >= 0.0  # intersects upper endcap
		return Ray([0.0,0.0,  cl.zmax], [0.0,0.0,1.0])
	elseif t == 0 && r.d[3] < 0.0  # intersects lower endcap
		return Ray([0.0,0.0, cl.zmin], [0.0, 0.0, -1.0])
	end
    
	P = ray(r, t)
    z = P[3]
    if z < cl.zmin || z > cl.zmax
        t, P = ray_intersection_with_cylinder_end_caps(r, cl, t)
		return Ray(P, r.d)
	else
		return Ray(P, r.d)
	end
end

# ╔═╡ 24cd55d7-6752-4e48-a3fa-7951cbff1b36
 rxp = ray_intersection_with_cylinder(rzp, cyl)

# ╔═╡ c3abda05-42ec-4759-bea4-cdf7ee92ed70
in_endcaps(cyl, rxp.e)

# ╔═╡ 366cf831-3558-44f2-bf09-6d9de74bf04d
 rxn = ray_intersection_with_cylinder(rzn, cyl)

# ╔═╡ 95a82d89-7b7f-4ce9-8f21-b527543c57f7
in_endcaps(cyl, rxn.e)

# ╔═╡ 21aeb736-8655-4747-bfa7-1ad9b5395fc3
 rxy = ray_intersection_with_cylinder(ry, cyl)

# ╔═╡ 5dcb797a-d7ef-463b-a629-336a8752c4ed
in_endcaps(cyl, rxy.e)

# ╔═╡ 92ef9d4b-ebf7-4247-a13e-aebe60cd79ad
cylinder_equation(cyl, rxy.e)

# ╔═╡ 80b17191-b750-4f31-9577-1f5214440429
in_barrel(cyl, rxy.e)

# ╔═╡ a143fe7e-e647-422e-8e0b-65099ff14d88
 rxz = ray_intersection_with_cylinder(rz, cyl)

# ╔═╡ 0dc2ca61-8ca9-43f2-ba71-87dcc2be08ec
in_endcaps(cyl, rxz.e)

# ╔═╡ 312154b9-aa06-4413-bcd7-7688fe168660
in_barrel(cyl, rxz.e)

# ╔═╡ 24cc705d-07dd-49bf-9542-49ce3b7c4855
"""
Transport gammas through the copper shield

"""
function transport_gammas_cs(c::Cylinder, df::DataFrame)
	nphotons = size(df)[1] # these are the points in df
    
    vex = zeros(nphotons)
	vey = zeros(nphotons)
	vez = zeros(nphotons)
	gex = zeros(nphotons)
	gey = zeros(nphotons)
	gez = zeros(nphotons)
	gdx = zeros(nphotons)
	gdy = zeros(nphotons)
	gdz = zeros(nphotons)

    for i in 1:nphotons
    	d = [df[i,"gdx"], df[i,"gdy"], df[i,"gdz"]] # direction
		p = [df[i,"gex"], df[i,"gey"], df[i,"gez"]] # position
        r = Ray(p,d)
        irr = ray_intersection_with_cylinder(r, c)
		vex[i] = p[1]
		vey[i] = p[2]
		vez[i] = p[3]
        gex[i] = irr.e[1]
		gey[i] = irr.e[2]
		gez[i] = irr.e[3]
		gdx[i] = irr.d[1]
		gdy[i] = irr.d[2]
		gdz[i] = irr.d[3]
	end
    DataFrame(vex=vex, vey=vey, vez=vez,
		      gex=gex, gey=gey, gez=gez,
		      gdx=gdx, gdy=gdy, gdz=gdz)
end

# ╔═╡ 14583522-a4d2-4f37-9c8b-5e2f54ace834
"""
Generate three standard normally distributed numbers and normalize the vector.

We have to be careful in the case that the vector has a norm close to zero, in which we must worry about floating point precision by dividing by a very small number. This is the reason for the while loop.
"""
function vectors_spherical(npoints::Integer, eps=1e-4; seed=123)

	Random.seed!(seed)
	mean = zeros(3)
	C = ScalMat(3, 1.0) 
	d = MvNormal(mean, C)

	vx = zeros(npoints)
	vy = zeros(npoints)
	vz = zeros(npoints)
	
	for i = 1:npoints
		v = zeros(3)
	    while norm(v) < eps
	        v = rand(d, 1)
	    end
    
    	v = v / norm(v)
		vx[i] = v[1]
		vy[i] = v[2]
		vz[i] = v[3]
	end
	 DataFrame(ID=1:npoints, vx=vx, vy=vy, vz=vz)
end

# ╔═╡ 934e0612-9284-41f5-a58c-03e002094366
"""
Generate three standard normally distributed numbers and normalize the vector.

We have to be careful in the case that the vector has a norm close to zero, in which we must worry about floating point precision by dividing by a very small number. This is the reason for the while loop.
"""
function vectors_spherical2(npoints::Integer, eps=1e-5)

	vx = zeros(npoints)
	vy = zeros(npoints)
	vz = zeros(npoints)

	for i = 1:npoints
		v = zeros(3)
	    while norm(v) < eps
	        x = randn()  # random standard normal
	        y = randn()
	        z = randn()
	        v = [x, y, z]
	    end
    
    	v = v / norm(v)
		vx[i] = v[1]
		vy[i] = v[2]
		vz[i] = v[3]
	end
	 DataFrame(ID=1:npoints, vx=vx, vy=vy, vz=vz)
end

# ╔═╡ 1f6b02e3-ad55-44f3-a417-39285058aee5
gdf = vectors_spherical2(1000)

# ╔═╡ c80fa6be-29d5-4296-982a-f28a7a5bba90
@rtransform gdf begin
	:norm = norm([:vx, :vy, :vz])
end

# ╔═╡ cef79ccd-b9ac-4a63-9729-bfbc23906c6e
GLMakie.scatter(gdf.vx, gdf.vy, gdf.vz)

# ╔═╡ f55a22e9-4d24-490d-a4f7-5b9ca0f5cda9
"""
gtransport is a short name for
`generate photons in point p inside cylinder and propagate to cylinder surface'

"""
function gtransport(c::Cylinder, p::Vector{Float64}, nphotons::Integer=10)

    rdf = vectors_spherical2(nphotons) #
    gex = zeros(nphotons)
	gey = zeros(nphotons)
	gez = zeros(nphotons)
	gdx = zeros(nphotons)
	gdy = zeros(nphotons)
	gdz = zeros(nphotons)

    for i in 1:nphotons
    	d = [rdf[i,"vx"], rdf[i,"vy"], rdf[i,"vz"]]
        r = Ray(p,d)
        irr = ray_intersection_with_cylinder(r, c)
        gex[i] = irr.e[1]
		gey[i] = irr.e[2]
		gez[i] = irr.e[3]
		gdx[i] = irr.d[1]
		gdy[i] = irr.d[2]
		gdz[i] = irr.d[3]
	end
    gdf = DataFrame(id=1:nphotons, gex=gex, gey=gey, gez=gez, 
		                           gdx=gdx, gdy=gdy, gdz=gdz)
	leftjoin(rdf, gdf, on="ID"=>"id", matchmissing=:equal)
	#rdf,gdf
end

# ╔═╡ 47ff0204-6b05-4635-89fa-23dac3a68500
gtdf = gtransport(stdsk, gp, np)

# ╔═╡ a0919f35-c839-4bf7-a368-ae7a165acc74
dfs =dropmissing(@chain gtdf begin                        
	   @rsubset in_barrel(stdsk, [:gex, :gey, :gez]) == true
	   @transform(_, :phi = atan.(:gey, :gex) .+ π, :z = :gez)
	   @select(_, :phi, :z)
	   #@rtransform(_, :gex = 50 + :gex)
	   #@select!  (:gex, :gey, :gez)
	   #transform(_, [:gex]  => xindx => :ix)
	   #transform(_, [:gey]  => xindx => :iy)
	end)

# ╔═╡ 952f6746-1cbb-4db0-9b5c-524a40caa787
gexmin = minimum(dfs.phi)

# ╔═╡ be1db55d-6c1f-4772-93cc-774702b7f5cb
geymin = minimum(dfs.z)

# ╔═╡ 606dca3d-d129-4ce1-86a6-ac9d05575c98
begin
	ffz = Figure()
	stephist(ffz[1, 1], dfs.phi, bins = 20)
	stephist(ffz[1, 2], dfs.z, bins = 20)
	ffz
end

# ╔═╡ ac6c096e-5085-436b-9cfb-5855d2312acb
 dfs2 = @chain dfs begin                        
		   transform(_, [:phi]  => pindx => :iphi)
		   transform(_, [:z]  => zindx => :iz)
	end

# ╔═╡ 223f2e63-71d1-44e7-9d30-925435251944
# group in (ix, iy), count the number of rows per group 
	# (which gives the number of photons per cell)
	# take the energy of the cell as the sum of the photons x pdf (SiPM)
	# take x,y as the barycenter of positions
	dfbx = unique(combine(groupby(dfs2, [:iphi, :iz])) do xypdf
				nb = nrow(xypdf) 
				DataFrame(
					      iphi =xypdf.iphi,
					      iz =xypdf.iz,
				          occ=nb/np)
    end)

# ╔═╡ da151c13-df1f-4d23-95de-93ecffbed3e3
stephist(dfbx.occ, bins = 20)

# ╔═╡ 610fdd44-266d-4d1f-8c20-d6ad6ddb6af3
"""
Transport gammas in dataframe df to the surface of the cylinder

"""
function transport_gammas(c::Cylinder, df::DataFrame)
	nphotons = size(df)[1] # these are the points generated in cylinder
    
    vex = zeros(nphotons)
	vey = zeros(nphotons)
	vez = zeros(nphotons)
	gex = zeros(nphotons)
	gey = zeros(nphotons)
	gez = zeros(nphotons)
	gdx = zeros(nphotons)
	gdy = zeros(nphotons)
	gdz = zeros(nphotons)

    for i in 1:nphotons
		rdf = vectors_spherical2(1) # throw just one random
    	d = [rdf[1,"vx"], rdf[1,"vy"], rdf[1,"vz"]] # direction
		p = [df[i,"gex"], df[i,"gey"], df[i,"gez"]]
        r = Ray(p,d)
        irr = ray_intersection_with_cylinder(r, c)
		vex[i] = p[1]
		vey[i] = p[2]
		vez[i] = p[3]
        gex[i] = irr.e[1]
		gey[i] = irr.e[2]
		gez[i] = irr.e[3]
		gdx[i] = irr.d[1]
		gdy[i] = irr.d[2]
		gdz[i] = irr.d[3]
	end
    DataFrame(vex=vex, vey=vey, vez=vez,
		      gex=gex, gey=gey, gez=gez,
		      gdx=gdx, gdy=gdy, gdz=gdz)
end

# ╔═╡ 96f53108-a469-4449-9563-0d4e193cebcc
"""
Generate points inside a cylinder.
"""
function points_in_cylinder(cyl::Cylinder, npoints::Integer; seed=12345)
	Random.seed!(seed)
	r = cyl.r * sqrt.(rand(Float64, npoints)) # random points along r
	theta = 2π * rand(Float64, npoints) # random points along theta
	tr = 2π * rand(Float64, npoints)
	vx =  r .* cos.(theta)
	vy =  r .* sin.(theta)
	#vx = cyl.r * rand(Uniform(-1,1), npoints)
	#vy = cyl.r * rand(Uniform(-1,1), npoints)
	vz = clength(cyl) * rand(Float64, npoints)
	DataFrame(ID=1:npoints, gex=vx, gey=vy, gez=vz)
end

# ╔═╡ d9ef36e2-56f6-4655-98f2-37c2699363b5
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

# ╔═╡ 32b53e87-034e-49db-842f-58f0c0002ed7
"""Normal to the cylinder barrel
Uses equation of cylynder: F(x,y,z) = x^2 + y^2 - r^2 = 0
then n = Grad(F)_P /Norm(Grad(F))_P
n = (2x, 2y, 0)/sqrt(4x^2 + 4y^2) = (x,y,0)/r (P)
"""
normal_to_barrel(c::Cylinder, P::Vector{Float64}) = [P[1], P[2], 0] ./ c.r



# ╔═╡ cb1d79b8-ada9-481e-8169-ab5976cb0759
"""
Give two points p0 and p1, thid function recturns;
mag: the magnitude of the vector defining the axis between p1 and p2
v: the unit vector in the direction of the axis defined by p0 and p1
n1, an n2: two unit vectors perpendicular to v and among themselves. 
"""
function cyl_unit_vectors(p0::Vector{Float64}, p1::Vector{Float64})
	#vector in direction of axis
	v = p1 - p0

	#find magnitude of vector
	mag = norm(v,2)

	#unit vector in direction of axis
	v = v / mag

	# choose (1,0,0) as second axis unless is first axis
	not_v = [1.0, 0.0, 0.0]
	if v == not_v
		not_v = [0.0, 1.0, 0.0]
	end

	#make vector perpendicular to v and not v
	n1 = cross(v, not_v)

	#normalize n1
	n1 = n1 /norm(n1,2)

	#make unit vector perpendicular to v and n1
	n2 = cross(v, n1)

	mag, v,  n1, n2
end

# ╔═╡ fdcde3b3-f80d-4852-9734-aa3b232e6b1d
"""
Parameterize the three surfaces of a cylinder
"""
function cylinder_surfaces_(c::Cylinder, np, nz, nr)

	mag, v,  n1, n2 = cyl_unit_vectors(c.p0, c.p1)
	
	#surface ranges over t from 0 to length of axis and 0 to 2*pi
	t = collect(range(0., mag, nz))
	theta = collect(range(0., 2π, np)) 
	rsample = collect(range(0., c.r, 2)) 

	#use meshgrid to make 2d arrays
	t, theta2 = jmgrid(t, theta)

	rsample,theta = jmgrid(rsample, theta)

	#generate coordinates for surface
	# "Tube"
	#X,Y,Z = [c.p0[i] .+ v[i] * t .+ c.r * sin.(theta2) * n1[i] .+ c.r * cos.(theta2) *  n2[i] for i in 1:3]
	
	X,Y,Z = [v[i] * t .+ c.r * sin.(theta2) * n1[i] .+ c.r * cos.(theta2) *  n2[i] for i in 1:3]
	
	# "Bottom"
	
	#X2, Y2, Z2 = [c.p0[i] .+ rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	X2, Y2, Z2 = [rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	
	# "Top"
	#X3, Y3, Z3 = [c.p0[i] .+ v[i]*mag .+ rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	X3, Y3, Z3 = [v[i]*mag .+ rsample .* sin.(theta) * n1[i] .+ rsample .* cos.(theta) * n2[i] for i in 1:3]
	
	t, theta2, theta, rsample, (X,Y,Z), (X2, Y2, Z2), (X3, Y3, Z3)
end


# ╔═╡ 5a74d977-a44e-4cba-b738-7bdd76056686
"""
Shows explicitly how to parameterize and draw the cylinder barrel
"""
function draw_cylinder_barrel(c::Cylinder, ptheta=100, pz=100)
	r = c.r
	h = clength(c)
	m, n =ptheta, pz
	u = range(0, 2π , length=n)
	v = range(0, h, length=m)
	
	us = ones(m)*u'
	vs = v*ones(n)'
	#Surface parameterization
	X = r*cos.(us)
	Y = r*sin.(us)
	Z = vs
	GLMakie.surface(X, Y, Z)
end

# ╔═╡ 53a6f956-1957-4915-9900-3af1f44ed2fa
draw_cylinder_barrel(cyl, 2, 100)

# ╔═╡ 252b84b2-aa7c-4d17-887b-241d0302e446
"""
Draw the three surfaces of a cylinder (barrel, bottom, top) using GLMakie
"""
function draw_cylinderx(cl::Cylinder, np=100, nz=10, 
	                    nr=10; viewa=[1,0,1], viewr=π,
                        showbarrel=true, showendcups=true)


	# Get the surfaces 
	ts,theta2s,thetas, rs, P1, P2, P3 = cylinder_surfaces_(cl, np, nz, nr)
	
	scene = Scene()
	cam3d!(scene)
	
	if showbarrel 
		sfs = GLMakie.surface!(scene, P1..., transparency = true, fxaa = true, overdraw=false, colormap=:reds)
	end
	if showendcups
		sfs2 = GLMakie.surface!(scene, P2..., transparency = true, fxaa = true, overdraw=false, colormap=:blues)
		sfs3 = GLMakie.surface!(scene, P3..., transparency = true, fxaa = true, overdraw=true, colormap=:greens)
	end
	
	#GLMakie.rotate!(scene, Vec3f(1, 0, 1), 1.5)
	GLMakie.rotate!(scene, Vec3f(viewa...), viewr)
	center!(scene)
	scene
end

# ╔═╡ 911d4179-909c-48f5-9a40-f23efe044498
draw_cylinderx(cyl, 200) #draws the three surfaces (cylinder is transparent)

# ╔═╡ 646c967d-145a-4f2a-9330-c4dbca8d78fc
draw_cylinderx(cyl, 200, showbarrel=true, showendcups=false,viewa=[1,0,0], viewr=2π)

# ╔═╡ 4a047576-ab00-4821-9383-1ba945f39285
draw_cylinderx(cyl, 200, showbarrel=false, showendcups=true)

# ╔═╡ 476bc4d8-7e22-4342-92ad-4896e47d437f
draw_cylinderx(stdsk, 200)

# ╔═╡ 72a529b6-01b2-44c4-b24c-139877cbb5d3
"""
Draw the three surfaces of a cylinder (barrel, bottom, top) using GLMakie
"""
function draw_cylinder2x(cl1::Cylinder, cl2::Cylinder, np=100, nz=10, 
	                    nr=10; viewa=[1,0,1], viewr=π,
                        showbarrel=true, showendcups=true)


	# Get the surfaces 
	_,_,_, _, P11, P21, P31 = cylinder_surfaces_(cl1, np, nz, nr)
	_,_,_, _, P12, P22, P32 = cylinder_surfaces_(cl2, np, nz, nr)
	scene = Scene()
	cam3d!(scene)
	
	if showbarrel 
		sfs = GLMakie.surface!(scene, P11..., transparency = true, fxaa = true, overdraw=false, colormap=:reds)
		sfs = GLMakie.surface!(scene, P12..., transparency = true, fxaa = true, overdraw=false, colormap=:blues)
	end
	if showendcups
		sfs2 = GLMakie.surface!(scene, P21..., transparency = true, fxaa = true, overdraw=false, colormap=:rainbow)
		sfs2 = GLMakie.surface!(scene, P21..., transparency = true, fxaa = true, overdraw=false, colormap=:greens)
		sfs3 = GLMakie.surface!(scene, P31..., transparency = true, fxaa = true, overdraw=true, colormap=:grays)
		sfs3 = GLMakie.surface!(scene, P32..., transparency = true, fxaa = true, overdraw=true, colormap=:bluesreds)
	end
	
	#GLMakie.rotate!(scene, Vec3f(1, 0, 1), 1.5)
	GLMakie.rotate!(scene, Vec3f(viewa...), viewr)
	center!(scene)
	scene
end

# ╔═╡ a41e2bc7-533c-40f7-b69f-94171750a44a
"""
Draw the three surfaces of a cylinder (barrel, bottom, top) using GLMakie
"""
function draw_rays_cylinderx(cl::Cylinder, gtdf, np=100, nz=10, 
	                    nr=10; viewa=[1,0,1], viewr=π,
                        showbarrel=true, showendcups=true, showrays=true)


	# Get the surfaces 
	ts,theta2s,thetas, rs, P1, P2, P3 = cylinder_surfaces_(cl, np, nz, nr)
	
	scene = Scene()
	cam3d!(scene)

	if showbarrel 
		sfs = GLMakie.surface!(scene, P1..., transparency = true, fxaa = true, overdraw=false, colormap=:reds)
	end
	if showendcups
		sfs2 = GLMakie.surface!(scene, P2..., transparency = true, fxaa = true, overdraw=false, colormap=:blues)
		sfs3 = GLMakie.surface!(scene, P3..., transparency = true, fxaa = true, overdraw=true, colormap=:greens)
	end
	if showrays
		GLMakie.scatter!(scene, gtdf.gex, gtdf.gey, gtdf.gez, color=:black, overdraw=false)
	end
	
	GLMakie.rotate!(scene, Vec3f(viewa...), viewr)
	center!(scene)
	scene
end

# ╔═╡ 09ec4f0d-bd08-4fee-83b4-40844a2b7236
draw_rays_cylinderx(stdsk, gtdf)

# ╔═╡ 7dbeff2f-2b28-4adf-b79a-d6b9b8c000d5
draw_rays_cylinderx(stdsk, gtdf, viewa=[1,0,0], viewr=2π)

# ╔═╡ 78059f1f-5056-4da5-94b3-f995413ebbe4


# ╔═╡ 2f45109c-70d4-4050-99d8-5c999fb94ec3


# ╔═╡ 73d06337-6871-4e10-8c80-fac2346f8db4
function lines_in_3D()
    #seed!(123)
    n = 10
    x, y, z = randn(n), randn(n), randn(n)
	println("x =",x)
    aspect=(1, 1, 1)
    perspectiveness=0.5
    # the figure
    fig = Figure(; resolution=(1200, 500))
    ax1 = Axis3(fig[1, 1]; aspect, perspectiveness)
    ax2 = Axis3(fig[1, 2]; aspect, perspectiveness)
    ax3 = Axis3(fig[1, 3]; aspect=:data, perspectiveness)
    lines!(ax1, x, y, z; color=1:n, linewidth=3)
    scatterlines!(ax2, x, y, z; markersize=15)
    hm = meshscatter!(ax3, x, y, z; markersize=0.2, color=1:n)
    lines!(ax3, x, y, z; color=1:n)
    #Colorbar(fig[2, 1], hm; label="values", height=15, vertical=false,
     #   flipaxis=false, ticksize=15, tickalign=1, width=Relative(3.55/4))
    fig
end

# ╔═╡ 56e32497-e335-4ca4-b1c8-3ec9ebcef24f
lines_in_3D()

# ╔═╡ Cell order:
# ╠═a298b03b-c7c4-4114-9984-d37e7df8d4c8
# ╠═5158389a-fca1-11ed-05a1-c7345c9bd354
# ╠═e24e5710-04b9-49c3-b3f6-2b136a1744eb
# ╠═e36751ff-a521-4529-95c8-6bbfc3314e66
# ╠═df3dde7a-9437-409f-bd51-db975f74b78c
# ╠═fc869138-1cf7-463b-aca5-729741f8b1ff
# ╠═f230deec-f4c9-429a-89f5-92d8c08e8472
# ╠═7371f833-a7a8-482c-925b-730d8378c22e
# ╠═2c803aef-9a52-4108-a252-e63a95929c67
# ╠═c7583c88-4cfc-49e1-b61e-c299801be450
# ╠═3c398d08-ef94-4202-9b87-62aa4fba07ed
# ╠═53a6f956-1957-4915-9900-3af1f44ed2fa
# ╠═ef077643-f5b1-4244-af7e-8111ede962a0
# ╠═81d62300-1a88-4c9c-b050-8acd20857867
# ╠═5e1c520f-fb11-4698-a981-5f3863c4788c
# ╠═911d4179-909c-48f5-9a40-f23efe044498
# ╠═646c967d-145a-4f2a-9330-c4dbca8d78fc
# ╠═4a047576-ab00-4821-9383-1ba945f39285
# ╠═c48c5f53-2632-4cfc-b07b-7fde9f1d8eb6
# ╠═de0fb8fd-c757-4769-a2da-d874add0a553
# ╠═50f7c64a-4fa6-4705-b2e9-896ef545029f
# ╠═d06eaf38-f8bb-4732-9b8f-92f2a1053176
# ╠═2940f9d9-f2c3-4772-8208-2375370b7e61
# ╠═a7378f23-9b56-4a5f-8c4b-b51545e69633
# ╠═c95dcf93-40ff-4db0-a323-29991f5135dd
# ╠═c7200fa5-31eb-458c-80ca-70f968edd9b0
# ╠═9e193e81-0147-46d2-943f-1526dfb26e50
# ╠═a89c0c98-ee2c-4103-8cd3-9be03a2b75ae
# ╠═da50ac8a-8381-490b-a000-b6ba8b392888
# ╠═aa445ed2-dc2d-4802-ab12-18505487ceef
# ╠═24cd55d7-6752-4e48-a3fa-7951cbff1b36
# ╠═6c2b1465-02ae-4768-a4f6-e8f461339b97
# ╠═c3abda05-42ec-4759-bea4-cdf7ee92ed70
# ╠═26d72ea7-31e5-4c55-af55-7666f7c47f1c
# ╠═b9c450a6-446e-4a5c-93ee-28a45cf96356
# ╠═5de835fa-284f-4dd0-ba1d-033d4c50da5a
# ╠═71203991-d5c4-439f-8dd1-fe1b4cf3b875
# ╠═366cf831-3558-44f2-bf09-6d9de74bf04d
# ╠═95a82d89-7b7f-4ce9-8f21-b527543c57f7
# ╠═59138131-46c1-4ac8-be0b-42f5ca2e30b8
# ╠═02a7fd08-440c-4e0c-93c8-d42792642d11
# ╠═40cff386-97de-4e6c-8d22-216302831893
# ╠═bca3de2e-74a6-423e-a660-3845c3e9657b
# ╠═7f4ab967-b5c2-46d9-90c4-ebbbdab0dbfe
# ╠═b78b7410-886f-47db-b8a3-b8a7feeec341
# ╠═c2d3e54d-3ad5-46b6-9d28-be8acf1086b2
# ╠═96af216d-c0e3-4825-a11f-3aaa04c42d0b
# ╠═21aeb736-8655-4747-bfa7-1ad9b5395fc3
# ╠═5dcb797a-d7ef-463b-a629-336a8752c4ed
# ╠═92ef9d4b-ebf7-4247-a13e-aebe60cd79ad
# ╠═80b17191-b750-4f31-9577-1f5214440429
# ╠═c97578f2-1d39-497a-9e79-31e041e77b7b
# ╠═4d17bf89-d01b-4f1b-9dc7-9dcd00ee1d8a
# ╠═d137c26c-bf75-4899-ac0a-eb8e70c792bb
# ╠═a143fe7e-e647-422e-8e0b-65099ff14d88
# ╠═0dc2ca61-8ca9-43f2-ba71-87dcc2be08ec
# ╠═312154b9-aa06-4413-bcd7-7688fe168660
# ╠═68e9ab79-21a4-48c1-82c3-45d8862e70cf
# ╠═1f6b02e3-ad55-44f3-a417-39285058aee5
# ╠═b4110420-ef1c-4e1b-868a-68bd9d943fea
# ╠═c80fa6be-29d5-4296-982a-f28a7a5bba90
# ╠═5e866260-a798-4c5f-b407-38d091d0ed07
# ╠═cef79ccd-b9ac-4a63-9729-bfbc23906c6e
# ╠═be8dd99e-60f5-4d75-959c-cffbd5bfa4b8
# ╠═cf2cce55-84fa-408a-b7ff-6748754eea34
# ╠═c7d4775a-f585-4852-a172-337655cb21ed
# ╠═476bc4d8-7e22-4342-92ad-4896e47d437f
# ╠═8a45b62d-cbfb-4f8c-9c61-fad0e3ca53da
# ╠═f99d7b04-2fc6-48c5-a846-601e2e270691
# ╠═128d5ed6-3a09-427e-b50d-159d0288d1dd
# ╠═47ff0204-6b05-4635-89fa-23dac3a68500
# ╠═09ec4f0d-bd08-4fee-83b4-40844a2b7236
# ╠═7dbeff2f-2b28-4adf-b79a-d6b9b8c000d5
# ╠═8baa74f5-557e-49cb-a6c2-a18be5e3c1d1
# ╠═63bedd95-3e1b-4786-8a14-91db66e503fd
# ╠═31628f89-3dcc-461a-81b2-582420805f91
# ╠═571fe2c8-4b2b-4436-b84b-f59bf9d32b13
# ╠═b6a5d9b3-60d2-44b9-ab35-cb298466b0d4
# ╠═a0919f35-c839-4bf7-a368-ae7a165acc74
# ╠═952f6746-1cbb-4db0-9b5c-524a40caa787
# ╠═be1db55d-6c1f-4772-93cc-774702b7f5cb
# ╠═35f50d8e-e56a-4e9a-8288-e15b133f065e
# ╠═606dca3d-d129-4ce1-86a6-ac9d05575c98
# ╠═ac6c096e-5085-436b-9cfb-5855d2312acb
# ╠═223f2e63-71d1-44e7-9d30-925435251944
# ╠═da151c13-df1f-4d23-95de-93ecffbed3e3
# ╠═5cca07ea-8740-4899-bca7-629bcdb8d46b
# ╠═5db4cd9d-71bf-4614-b6b6-922f3dcbdbf1
# ╠═953f8162-fb20-4062-b185-3271340a135f
# ╠═ce9da661-5f2f-4e5d-825e-5d029c274066
# ╠═24b1b8ef-10ee-4d96-8121-827234233873
# ╠═6458b9a6-e0b1-40cf-87bd-3e3f01990525
# ╠═478e01bc-d764-4af7-a6f4-661cdba4820d
# ╠═efe04df4-cb56-4916-96f7-296dc37492ac
# ╠═2f01a131-d184-4ea8-97c1-89309716383f
# ╠═053952a1-0558-447c-a3f5-ae30c2d2f7be
# ╠═991a026a-bfe2-4bc8-a200-31dd25a881d1
# ╠═2facda59-dc5a-4951-8d6f-1d121690aad2
# ╠═539faa01-299b-4c5b-8113-fc21501dd84d
# ╠═269e4b09-3505-40fe-8815-2bc847d02a99
# ╠═f55a22e9-4d24-490d-a4f7-5b9ca0f5cda9
# ╠═610fdd44-266d-4d1f-8c20-d6ad6ddb6af3
# ╠═24cc705d-07dd-49bf-9542-49ce3b7c4855
# ╠═c26c8511-afae-4447-9458-d1d21f565911
# ╠═deefebf9-07a3-4cce-b2d5-4b42644b5db7
# ╠═6f6cd50f-d4fd-48d8-8b53-7edf9112c45b
# ╠═e9d59365-ea79-4e97-b908-cbd528dc0804
# ╠═c66ca1b8-1ae8-4b53-b828-88cac5d8eae9
# ╠═da14195e-414e-40ab-aca1-83d2d68b9172
# ╠═0f55b2c3-db7b-4aad-bd8d-c678a5db0b02
# ╠═f41fbfdf-4cee-4b82-b087-e224c70b7f6b
# ╠═8da08d4a-f0cf-4b36-b634-b9f8716a3184
# ╠═7f8c81be-0f2d-4214-b27d-a7470534d867
# ╠═e1f26c48-0c0c-480c-9574-96b29b65ccf7
# ╠═372b5f7a-2919-4c73-8ce4-3ebf8ebb21b9
# ╠═834027c9-0a39-450f-9996-a1ec56f080a9
# ╠═9022e8da-fe82-4631-8e83-549dad785f9a
# ╠═8569dc46-642f-4709-9a78-0292f5aa1613
# ╠═1061c7ff-fbf4-4dd4-bfc3-fb2c38022d87
# ╠═0abeba90-ba75-4231-9dd8-d35400964c60
# ╠═464a3cf9-c030-4270-bae7-bf4339ed3fc8
# ╠═58491933-49ec-4b70-ab1f-2301657ed9e8
# ╠═4d53eaca-0cad-4f34-8205-216e9c0a16d9
# ╠═989cec01-548a-40e8-9a82-149dbfa7d367
# ╠═c6f63bc1-6d4d-49b3-988f-ccaaaae562e1
# ╠═c40ad59c-7687-4bba-b923-fa498295821c
# ╠═8547b5f7-6b8a-4a9b-b0bf-5b1d8b08c072
# ╠═14583522-a4d2-4f37-9c8b-5e2f54ace834
# ╠═934e0612-9284-41f5-a58c-03e002094366
# ╠═96f53108-a469-4449-9563-0d4e193cebcc
# ╠═d9ef36e2-56f6-4655-98f2-37c2699363b5
# ╠═32b53e87-034e-49db-842f-58f0c0002ed7
# ╠═cb1d79b8-ada9-481e-8169-ab5976cb0759
# ╠═fdcde3b3-f80d-4852-9734-aa3b232e6b1d
# ╠═5a74d977-a44e-4cba-b738-7bdd76056686
# ╠═252b84b2-aa7c-4d17-887b-241d0302e446
# ╠═72a529b6-01b2-44c4-b24c-139877cbb5d3
# ╠═a41e2bc7-533c-40f7-b69f-94171750a44a
# ╠═78059f1f-5056-4da5-94b3-f995413ebbe4
# ╠═2f45109c-70d4-4050-99d8-5c999fb94ec3
# ╠═73d06337-6871-4e10-8c80-fac2346f8db4
# ╠═56e32497-e335-4ca4-b1c8-3ec9ebcef24f
