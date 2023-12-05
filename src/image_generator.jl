using StaticArrays
using LinearAlgebra
using Printf
using Images
using Random
using DataFrames
using DataFramesMeta
using Logging
using CSV

"""
Defines a ray (e.g, a segment that propagates in straight line)
r = e + t * d, 
where e in the initial point and d is the direction vector
"""
struct Ray
    e::SVector{3,Float32}
    d::SVector{3,Float32}
end


"""
Propagates ray r along the distance defined by t
"""
ray(r::Ray, t::Float32) =r.e + t * r.d


"""
Defines a Plane

A x + B y + C  Z + D=0

where norm([A,B,C]) = 1
"""
struct Plane
    P::SVector{3,Float32} # [A, B, C]
    D::Float32
end

"""
Computes the intersection of a ray with a plane. 
Returns large if the ray does not intersect plane, t =v0/vd otherwise

"""
function ray_plane(r::Ray, p::Plane, eps=1e-7, large=1e+7)
	vd = dot(r.d, p.P)
	v0 = -(dot(r.e, p.P) + p.D)
	if abs(vd) < eps
		return large
	else
		return v0/vd
	end
end


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


"""
Propagate a ray to the rcub. It propagates to all surfaces
and returns the one with the shortest distance
"""
function propagate_ray_rcub(mcr::Rcub, rz::Ray, large=1.0f+7)
	dst = zeros(Float32, 6)
	dst[1] = ray_plane(rz, mcr.X0)
	dst[2] = ray_plane(rz, mcr.X1)
	dst[3] = ray_plane(rz, mcr.Y0)
	dst[4] = ray_plane(rz, mcr.Y1)
	dst[5] = ray_plane(rz, mcr.Z0)
    dst[6] = ray_plane(rz, mcr.Z1)

	I = []
	for i in 1:6
		@debug "i" i
		@debug "dst[i]" dst[i]
		if dst[i] > 0 && dst[i] < large
			append!(I,i)
		end
	end
	@debug "I" I
	xmin = large
	for ii in I
		if dst[ii] < xmin
			xmin = dst[ii]
		end
	end
	xmin
end




"""
Return dimensions of cuboid

"""
function dimensions(r::Rcub)
	(x0=-r.X0.D, x1=-r.X1.D, y0=-r.Y0.D, y1=-r. Y1.D,
     z0=-r.Z0.D, z1=-r.Z1.D)
end


"""
Create an Rcub from its three dimensiones

"""
function create_rcub(dx::Float32,dy::Float32,dz::Float32)
	X0 = Plane(SVector{3,Float32}(1,0,0), -0.0)
	X1 = Plane(SVector{3,Float32}(1,0,0), -dx)
	Y0 = Plane(SVector{3,Float32}(0,1,0), -0.0)
	Y1 = Plane(SVector{3,Float32}(0,1,0), -dy)
	Z0 = Plane(SVector{3,Float32}(0,0,1), -0.0)
	Z1 = Plane(SVector{3,Float32}(0,0,1), -dz)
	Rcub(X0,X1,Y0,Y1,Z0,Z1)
end


"""
Generate random points inside the cuboid

"""
function points_rcub(mcr::Rcub, npoints::Integer)
	dm = dimensions(mcr)

	vx = rand(dm.x0:0.1:dm.x1, npoints)
	vy = rand(dm.y0:0.1:dm.y1, npoints)
	vz = rand(dm.z0:0.1:dm.z1,npoints)
	DataFrame(ID=1:npoints, vx=Float32.(vx), vy=Float32.(vy), vz=Float32.(vz))
end


"""
Generate three standard normally distributed numbers 
and normalize the vector.

We have to be careful in the case that the vector 
has a norm close to zero, in which we must worry about 
floating point precision by dividing by a very small number. 
This is the reason for the while loop.

"""
function vectors_spherical2(npoints::Integer, eps=1e-5)
    
	vx = zeros(Float32, npoints)
	vy = zeros(Float32, npoints)
	vz = zeros(Float32, npoints)

	for i = 1:npoints
		v = vectors_spherical2(eps)
		vx[i] = v[1]
		vy[i] = v[2]
		vz[i] = v[3]
	end
	 DataFrame(ID=1:npoints, vx=vx, vy=vy, vz=vz)
end

"""
Overload the function with the version in which we generate a single vector
"""
function vectors_spherical2(eps=1e-5)
	v = zeros(3)
	while norm(v) < eps
		x = randn(Float32)  # random standard normal
		y = randn(Float32)
		z = randn(Float32)
		v = [x, y, z]
	end
    
    v / norm(v)
		
end


"""
Input:
 - A named tuple specifying the dimensions of the crystal
 - The number of optical photons to be generated
 - The transverse size of the SiPM in mm
 - The pde of the SiPM
 - rpos = true means that the point of generation is chosen randomly in the crystal
 - if rpos = false, the point of generation is xg, yg, zg 
 - returns an image of the Z detection plane as seen by SiPMs
"""
function ggun(dm::NamedTuple; ng::Integer=1000, xsipm::Float32=3.0, pde::Float32=1.0,
              rpos=true, xg::Float32=0.0f0, yg::Float32=0.0f0, zg::Float32=0.0f0,
              eps=1e-5)

	# Generates (x,y,z) inside crystal
	function getxyz(dm::NamedTuple)
		rand(dm.x0:0.1:dm.x1), rand(dm.y0:0.1:dm.y1), rand(dm.z0:0.1:dm.z1)
	end

	# copies extrapolation poin into p (exclusing outlayers)
	function fill_point!(P::SVector{3,Float32}, px::Vector{Float32}, 
	                     py::Vector{Float32},   pz::Vector{Float32}, 
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

	# Transforms a vector of values between 0.0 and dx (dy) in 
	# a vector of indices
	function indxy(xsipm::Float32)
		function indx(x::Vector{Float32}) 
			xx = Int.(floor.(x  ./xsipm ) )
			replace(xx, 0 =>1) 
		end
		indx
	end
	
	@debug "Crystal with dimensions" dm

    mcr = create_rcub(dm.x1, dm.y1, dm.z1)

    @debug "create rcub " mcr

	if rpos
		@debug "generate random position inside crystal"
		x,y,z = getxyz(dm) # generate position
	else
		@debug "Takes specified position"
		x, y, z = xg, yg, zg
	end
		@debug "Position" x y z

	@debug "Shooting optical photons =" ng
	
	#vdf =vectors_spherical2(ng) # direction vectors
	vdf =vectors_spherical2(ng, eps)
	typeof(vdf)
	#propagate gammas

	px = zeros(Float32, ng)
	py = zeros(Float32, ng)
	pz = zeros(Float32, ng)
	for (i, _) in enumerate(eachrow(vdf)) 
    	
		ry = Ray(SVector{3,Float32}(x,y,z), 
			     SVector{3,Float32}(vdf[i,"vx"], vdf[i,"vy"],vdf[i,"vz"]))	
		t = propagate_ray_rcub(mcr,ry)
		fill_point!(ray(ry,t), px, py, pz, i, dm)		
	end
	
	nbx = Int(floor(dm.x1 /xsipm))
	nby = Int(floor(dm.y1 /xsipm))
	@debug "Preparing an image of size" nbx nby	

	xindx = indxy(xsipm)  # generates indexes for SiPM of dimension xsipm
	
	rdf = DataFrame(ID=1:ng, vx=px, vy=py, vz=pz) #initial DF

	# Keep only photons reaching dz
	dfs = @chain rdf begin                        
		   @rsubset :vz == dm.z1 
		   @select!  $(Not(:vz))
		   transform(_, [:vx]  => xindx => :ix)
		   transform(_, [:vy]  => xindx => :iy)
	end

	# group in (ix, iy), count the number of rows per group 
	# (which gives the number of photons per cell)
	# take the energy of the cell as the sum of the photons x pdf (SiPM)
	# take x,y as the barycenter of positions
	dfbx = unique(combine(groupby(dfs, [:ix, :iy])) do xypdf
				nb = nrow(xypdf) 
				xbox = sum(xypdf.vx) / nb
				ybox = sum(xypdf.vy) / nb
				DataFrame(
					      ix =xypdf.ix,
					      iy =xypdf.iy,
				          ebox = nb * pde, xbox=xbox, ybox=ybox)
    end)

	exy = zeros(Float32, nbx,nby) #image
	#xr = zeros(nbx)
	#yr = zeros(nby)

	for r in eachrow(dfbx) 
		exy[r.ix,r.iy] = r.ebox
		#xr[r.ix] = r.xbox
		#yr[r.iy] = r.ybox
	end
	#DataFrame(eg=Int32(ng), xg=Float32(x), yg=Float32(y), zg=Float32(z)), exy ./ maximum(dfbx.ebox)  #gray scale, normalize to one
    DataFrame(ng=Int32(ng), xg=Float32(x), yg=Float32(y), zg=Float32(z)), 
    exy ./ maximum(dfbx.ebox)
end


"""
Uses ggun to shoot a number of gammas, propagate them to the SiPM
plane (Z1 plane) and create images and labels in specified directories.

"""
function imagerun(img0::Integer, imgl::Integer, dm::NamedTuple; 
                  ipath::String, lpath::String,
                  ng::Integer=50000, xsipm::Float32=6.0f0, pde::Float32=1.0f0,
                  rflbl = "lbl", rfimg = "img", rpos=true, seed=12345)

    function getlbls(i::Integer)
    flbl =string(rflbl,string(i))
    fimg =string(rfimg,string(i))
    plbl = joinpath(lpath, flbl)
    pimg = joinpath(ipath, fimg)
    pngimg = string(pimg,".png")
    csvlbl = string(plbl,".csv")
    csvlbl, pngimg
    end

    Random.seed!(seed)
    @info "saving images in directory $(ipath)" 
    @info "saving lbls in directory $(lpath)" 

    for i in img0:imgl
        csvlbl, pngimg = getlbls(i)
        @debug "csvlbl = $(csvlbl) pngimg = $(pngimg)"

        lbldf, imgx = ggun(dm; ng, xsipm, pde, rpos)
        @debug "label =>" lbldf
        @debug "image =>" Float32.(imgx)

        #imgmd = ImageMeta(imgx, xg=lbldf.xg[1], yg=lbldf.yg[1], 
        #                  zg=lbldf.zg[1], ng=lbldf.ng[1])

        @debug "saving image $i" 
        save(pngimg, imgx)

        @debug "saving labels $i" 
        CSV.write(csvlbl, lbldf)
    end
end
