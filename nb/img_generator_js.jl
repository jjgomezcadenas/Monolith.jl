
using Monolith
using Test
using StaticArrays
using Images
using CSV
using Glob
using DataFrames

import Monolith: Plane, Ray, Rcub, dimensions, create_rcub,
       ray_plane, propagate_ray_rcub,
       points_rcub, vectors_spherical2,
       ggun, imagerun, 
       load_images


# Import necessary functions
using FilePathsBase

function remove_old_file(dirpath)

	# List all files in the directory
	files = readdir(dirpath)

	# Optionally filter out only files, if you don't want to remove subdirectories
	#files = [file for file in files if isfile(joinpath(dir_path, file))]

	# Remove each file
	for file in files
		file_path = joinpath(dirpath, file)
		rm(file_path; force=true)
	end
end


# Data paths
rpath = "/Users/jjgomezcadenas/Data/monolith/testrun0"
ipath =joinpath(rpath,"imgs")
lpath =joinpath(rpath,"lbls")
	   
	   
# Dimensions of crystal
th = 18.6f0 * 2.0f0
dd = 50.0f0

# Size of SiPM
xsipm = 6.0f0

# PDE of SiPM
pde = 1.0f0

# Detector definition
mcr = create_rcub(dd,dd,th)
dm = dimensions(mcr)


nbx = Int(floor(dm.x1 /xsipm))
nby = Int(floor(dm.y1 /xsipm))
nbox = nbx * nby

@info "number of of boxes along x = ", nbx
@info "number of of boxes along y = ", nby
@info "number of of boxes         = ", nbox
	


lbldf, imgx = ggun(dm, ng=50000, xsipm=Float32(xsipm), 
                                       pde=Float32(pde), rpos=true)

#@info typeof(imgx)
@info "label data " lbldf


#tags = Dict{String,String}("xs"=>string(x), "y"=>string(y),
#"z"=>string(z))
#imgmd = ImageMeta(imgx, xg=x, yg=y,zg=z,eg=ng)

#save("test_img.png", imgmd)
#imgx2 = load("test_img.png")


#@info "array data " arraydata(imgmd)
#@info "meta data  " properties(imgmd)

# @info "writing file to disk"
# write_tags("test_img.jpg"; img=imgx, 
#                            tags=tags)

# @info "read back file"

# imgr = load("test_img.jpg")
# itags = read_tags("test_img.jpg")

# @info "recover image" imgr
# @info "recovered tags" itags
#imgmdr = ImageMeta(imgr, itags)

#@info "recovered array data " arraydata(imgmdr)
#@info "recovered meta data  " properties(imgmdr)

il = 1000
genimg = true 

rmfiles = true

if rmfiles
	@info "removing old files in $ipath"
	remove_old_file(ipath)

	@info "removing old files in $lpath"
	remove_old_file(lpath)
end

if genimg
	@info "Generating $il images "
	imagerun(1, il, dm; ng=500000, ipath=ipath, lpath=lpath,
	                    xsipm=Float32(xsipm), pde=Float32(pde), 
						rpos=true, seed=12345)

end
@info "Load example Image"

imgs = glob("*.png",ipath)
lbls = glob("*",lpath)
 
@info "reading test image $(imgs[1])" 
rimg = load(imgs[1])

@info "reading test label $(lbls[1])" 
rlbl =DataFrame(CSV.File(lbls[1]))

@info "image " Float32.(rimg)
@info "labels  " rlbl
# @info Float32.(rimg)

# lbls = glob("*",lpath)
# @info Arrow.Table(lbls[1])

