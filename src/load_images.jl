using DataLoaders
import DataLoaders.LearnBase: getobs, nobs
using Images
using Glob


function load_images_npy(ipath::String; train=1:2, norm=true)
	fnames = glob("*.npy",ipath)
	imgtrain =[]
	for it in train
		imgs = permutedims(load(fnames[it]), (3,2,1)) 
		imgx = zeros(Float32, size(imgs)...) 
		if norm
			for i in 1:size(imgs)[3]
				imgx[:,:,i] = Float32.(imgs[:,:,i]/maximum(imgs[:,:,i]))
			end
		end
		push!(imgtrain, imgx)
	end
	imgtrain
end