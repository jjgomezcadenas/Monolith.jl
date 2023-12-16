using Monolith
using Images
using ImageView
using CSV
using Glob
using DataFrames
import Monolith: load_images_npy

# Data paths
rpath = "/Users/jjgomezcadenas/Data/monolith/testrun1"
ipath =joinpath(rpath,"imgs")
lpath =joinpath(rpath,"lbls")
	   
imtrain = load_images_npy(ipath; train=1:2, norm=true)
imt1 = imtrain[1]
im1 = Gray.(imt1[:,:,3])
imshow(im1)


#gui = imshow_gui((8, 8), (2, 1))  # 2 columns, 1 row of images (each initially 300×300)
#canvases = gui["canvas"]
#imshow(canvases[1,1], Gray.(imt1[:,:,1]))
#show(gui["window"])
#imshow(canvases[1,2], Gray.(imt1[:,:,2]))