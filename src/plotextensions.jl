# This file needs a big overhaul in order to work properly!

import Gadfly

import Gadfly.cm
import Gadfly.pt


Gadfly.set_default_plot_size(15cm,15cm)

flipit(A) = flipud(A')

function plot(d::AbstractDomain{2}; n = 100)
	mask = in((n,n), d)
	Gadfly.spy(flipit(mask), Gadfly.Guide.xlabel("x"), Gadfly.Guide.ylabel("y"))
end

#function plot_extension{T}(f::ExpFun{1,T})
#	layer = lastlayer(f)
#	z = eval_extension(layer, layer.param.l)
#	grid = layer.context.gridL
#	Gadfly.plot(x=range(grid,1), y=z, Gadfly.Geom.line, Gadfly.Theme(line_width=4pt,minor_label_font_size=16pt,major_label_font_size=16pt))
#end

function plot{T}(f::ExpFun{1,T}; n=100)
	z,grid = evaluate(f, (n,))
	Gadfly.plot(x=range(grid,1), y=z, Gadfly.Geom.line, Gadfly.Theme(line_width=4pt,minor_label_font_size=16pt,major_label_font_size=16pt))
end

# function plot{T}(f::ExpFun{2,T}; n=100)
# 	z,grid = evaluate(f, (n,n,))
# 	Gadfly.spy(flipit(z), Gadfly.Guide.xlabel("x"), Gadfly.Guide.ylabel("y"))
# end


# import PyPlot

# function plot(d::AbstractDomain{2}; n = 100)
# 	mask = ExpFun.eval(d, (n,n))
# 	PyPlot.matshow(mask)
# 	PyPlot.tick_params(labelbottom="off",labeltop="off",labelleft="off",labelright="off")
# end

# function plot(f::ExpFun{2}; n = 100)
# 	z,grid = eval(f, (n,n))
# 	PyPlot.surf(range(grid,1), range(grid,2), z)
# end

# using PyCall

# @pyimport sip
# sip.setapi("QString", 2)
# sip.setapi("QVariant", 2)
# pygui_start(:qt)
# @pyimport mayavi.mlab as mlab


# function plot(f::ExpFun{2}; n = 100)
# 	z,grid = evaluate(f, (n,n))
# 	mlab.surf(range(grid,1), range(grid,2), z, representation="wireframe", warp_scale="auto")
# 	mlab.axes(x_axis_visibility=true, xlabel="x", ylabel="y")
# end

# function plot_extension(f::ExpFun{2}; n = 100)
# 	layer = lastlayer(f)
# 	z,grid = eval_extension(layer, (n,n))
# 	mlab.surf(range(grid,1), range(grid,2), z, representation="wireframe", warp_scale="auto")
# 	mlab.axes(x_axis_visibility=true, xlabel="x", ylabel="y")
# end

# function plot_accuracy{T}(f::ExpFun{2,T}, f_orig; n = 100)
# 	layer = lastlayer(f)
# 	z,grid = eval_extension(layer, (n,n))
# 	z2 = T[ f_orig(grid[i,j]...) for i=1:size(grid,1),j=1:size(grid,2) ]
# 	Gadfly.spy(flipit(log10(abs((z-z2)))), Gadfly.Guide.xlabel("x"), Gadfly.Guide.ylabel("y"))
# end

# function plot(d::AbstractDomain{3}; n=40)
# 	mask = ExpFun.in((n,n,n), d)
# 	mlab.contour3d(int(mask))
# end

# function plotslice(f::ExpFun{3}; n=100)
# 	z,grid = evaluate(f, (n,n,n))
# 	mlab.pipeline[:image_plane_widget](mlab.pipeline[:scalar_field](z), plane_orientation="x_axes", slice_index=10)
# 	mlab.pipeline[:image_plane_widget](mlab.pipeline[:scalar_field](z), plane_orientation="y_axes", slice_index=10)
# 	mlab.outline()
# end

# function plot(f::ExpFun{3}; n = 40)
# 	z,grid = call(f, (n,n,n))
# 	mlab.pipeline[:volume](mlab.pipeline[:scalar_field](z))
# end








