using DataLoaders
import DataLoaders.LearnBase: getobs, nobs
using Images
using Glob


#To use `DataLoaders`' multi-threading, you need to start Julia with multiple
# threads. Check the number of available threads with `Threads.nthreads()`.

struct ImageDataset
	files::Vector{String}
end


ImageDataset(folder::String) = ImageDataset(glob("*.png",folder))
	
nobs(data::ImageDataset) = length(data.files)
getobs(data::ImageDataset, i::Int) = Images.load(data.files[i])

function load_images(ipath::String)
	data = ImageDataset(ipath)
	for images in DataLoader(data, 1, collate = true)
	    @info images
	end
end