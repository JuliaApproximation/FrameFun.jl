# plots.jl

import BasisFunctions: plot_error, plot_expansion

# One-dimensional plot, just the domain
function plot(f::FrameFun{1}; n=201, color="blue")
    G = grid(resize(basis(f),n))
    G = MaskedGrid(G,domain(f))
    x = convert(Array{Float64},apply(x->x,eltype(f),G))
    data = convert(Array{Float64},real(f(G)))
    Main.PyPlot.plot(x,data, color=color)
    Main.PyPlot.title("Extension (domain)")
end

## # One-dimensional plot, including extension
## function plot_full{1}(f::Fun{1})

## end
function plot_expansion(f::FrameFun{1}; n=201, repeats=0, color="blue", alpha=1.0)
    G = grid(resize(basis(f),n))
    data = convert(Array{Float64},real(f(G)))
    x = convert(Array{Float64},apply(x->x,eltype(f),G))
    for i=-repeats:repeats
        Main.PyPlot.plot(x+i*(right(G)-left(G)),data,linestyle="dashed", color=color, alpha=alpha)
    end
    Main.PyPlot.plot(x,data, color=color, alpha=alpha)
    Main.PyPlot.title("Extension (Full)")
end

function plot_error(f::FrameFun{1}, g::Function; n=201, repeats = 0, color="blue")
    G = grid(resize(basis(f),n))
    data = real(f(G))
    x = convert(Array{Float64},apply(x->x,eltype(f),G))
    plotdata=convert(Array{Float64},abs(apply(g,eltype(f),G)-data))
    for i=-repeats:repeats
        Main.PyPlot.semilogy(x+i*(right(G)-left(G)),plotdata,linestyle="dashed", color=color)
    end
    Main.PyPlot.semilogy(x,plotdata,color=color)
    Main.PyPlot.ylim([min(minimum(log10(plotdata)),-16),1])
    Main.PyPlot.title("Absolute Error")
end

function plot_samples(f::FrameFun{1}; gamma=2)
    grid, fbasis2 = oversampled_grid(domain(f), basis(f), gamma)
    x = [grid[i] for i in eachindex(grid)]
    data = convert(Array{Float64},real(f(grid)))
    Main.PyPlot.stem(x,data)
    Main.PyPlot.title("samples")
end


## # Maybe place this in funs.jl?
## function call(f::FrameFun, g::AbstractGrid)
##     result = Array(eltype(f), size(g))
##     call!(f, result, g)
##     result
## end

## function call!{N}(f::FrameFun, result::AbstractArray, g::AbstractGrid{N})
##     x = Array(eltype(f), N)
##     for i in eachindex(g)
##         getindex!(x, g, i)
##         result[i] = call(f, x...)
##     end
## end

## function call!(f::FrameFun, result::AbstractArray, g::AbstractGrid1d)
##     for i in eachindex(g)
##         result[i] = call(f, g[i])
##     end
## end

## function call!(f::FrameFun, result::AbstractArray, x::AbstractArray)
##     @assert size(result) == size(x)
##     for i = 1:length(x)
##         result[i] = call(f, x[i])
##     end
## end


function apply(f::Function, return_type, g::AbstractGrid)
    result = Array(return_type, size(g))
    call!(f, result, g)
    result
end

function call!{N}(f::Function, result::AbstractArray, g::AbstractGrid{N})
    for i in eachindex(g)
        result[i] = f(getindex(g, i)...)
    end
end

function plot_domain(d::AbstractDomain{2}; n=1000)
    B = boundingbox(d)
    grid = equispaced_aspect_grid(B,n)
    Z = evalgrid(grid, d)
    Main.PyPlot.imshow(Z',interpolation="bicubic",cmap="Blues",extent=(left(B)[1], right(B)[1], left(B)[2], right(B)[2]),aspect="equal",origin="lower")
end


function plot(f::FrameFun{2};n=1000)
    B = boundingbox(domain(set(expansion(f))))
    Tgrid = equispaced_grid(B,n)
    Mgrid=MaskedGrid(Tgrid, domain(set(expansion(f))))
    data = convert(Array{Float64},real(expansion(f)(Mgrid)))
    x=[Mgrid[i][1] for i = 1:length(Mgrid)]
    y=[Mgrid[i][2] for i = 1:length(Mgrid)]
    Main.PyPlot.plot_trisurf(x,y,data)
end

function plot_image(f::FrameFun{2};n=300,unscaled=false, border=true, colorbar=true)
    d =domain(set(expansion(f)))
    Tgrid = grid(resize(basis(f), (n,n)))
    Mgrid=MaskedGrid(Tgrid, domain(f))
    Z = evalgrid(Tgrid, d)
    data = convert(Array{Float64},real(expansion(f)(Mgrid)))
    vmin = minimum(data)
    vmax = maximum(data)
    data = real(expansion(f)(Tgrid))
    if unscaled
        vmin = minimum(data)
        vmax = maximum(data)
    end
    X = convert(Array{Float64},real(apply((x,y)->x,eltype(f),Tgrid)))
    Y = convert(Array{Float64},real(apply((x,y)->y,eltype(f),Tgrid)))
    Main.PyPlot.pcolormesh(X,Y,data,vmin=vmin, vmax=vmax, shading = "gouraud")
    bound = boundary(Tgrid,d)
    border && plot_grid(bound)
    Main.PyPlot.axis("scaled")
    Main.PyPlot.xlim([left(Tgrid)[1], right(Tgrid)[1]])
    Main.PyPlot.ylim([left(Tgrid)[2], right(Tgrid)[2]])
    colorbar && Main.PyPlot.colorbar()
end

function plot_image(f::FrameFun{2}, g::Function; n=300, border=false)
    d =domain(set(expansion(f)))
    Tgrid = grid(resize(basis(f), (n,n)))
    Mgrid=MaskedGrid(Tgrid, domain(f))
    Z = evalgrid(Tgrid, d)
    data = real(apply(g,eltype(f),Mgrid))
    vmin = minimum(data)
    vmax = maximum(data)
    data = real(apply(g,eltype(f),Tgrid))
    X = convert(Array{Float64},real(apply((x,y)->x,eltype(f),Tgrid)))
    Y = convert(Array{Float64},real(apply((x,y)->y,eltype(f),Tgrid)))
    Main.PyPlot.pcolormesh(X,Y,data,vmin=vmin, vmax=vmax, shading = "gouraud")
    bound = boundary(Tgrid,d)
    border && plot_grid(bound)
    Main.PyPlot.axis("scaled")
    Main.PyPlot.xlim([left(Tgrid)[1], right(Tgrid)[1]])
    Main.PyPlot.ylim([left(Tgrid)[2], right(Tgrid)[2]])
    Main.PyPlot.colorbar()
end

function plot_error(f::FrameFun{2},g::Function;n=300,border=false)
    d =domain(set(expansion(f)))
    Tgrid = grid(resize(basis(f), (n,n)))
    Mgrid=MaskedGrid(Tgrid, domain(f))
    Z = evalgrid(Tgrid, d)
    data = real(expansion(f)(Mgrid))
    vmin = minimum(data)
    vmax = maximum(data)
    data = log10(abs(expansion(f)(Tgrid)-apply(g,eltype(f),Tgrid)))
    X = convert(Array{Float64},real(apply((x,y)->x,eltype(f),Tgrid)))
    Y = convert(Array{Float64},real(apply((x,y)->y,eltype(f),Tgrid)))
    Main.PyPlot.pcolormesh(X,Y,data,vmin=-16.0, vmax=1.0, shading = "gouraud")
    bound = boundary(Tgrid,d)
    border && plot_grid(bound)
    Main.PyPlot.axis("scaled")
    Main.PyPlot.xlim([left(Tgrid)[1], right(Tgrid)[1]])
    Main.PyPlot.ylim([left(Tgrid)[2], right(Tgrid)[2]])
    Main.PyPlot.colorbar()
    Main.PyPlot.title("log10 of absolute error")
end

function plot_grid(grid::AbstractGrid2d)
    x=[grid[i][1] for i = 1:length(grid)]
    y=[grid[i][2] for i = 1:length(grid)]
    Main.PyPlot.plot(x,y,linestyle="none",marker="o",color="black",mec="black",markersize=1.5)
    Main.PyPlot.axis("equal")
end

function plot_grid(grid::AbstractGrid3d)
    x=[grid[i][1] for i = 1:length(grid)]
    y=[grid[i][2] for i = 1:length(grid)]
    z=[grid[i][3] for i = 1:length(grid)]
    Main.PyPlot.plot3D(x,y,z,linestyle="none",marker="o",color="blue")
    Main.PyPlot.axis("equal")
end

function plot_grid{TG}(grid::TensorProductGrid{TG,2})
    dom = Cube(left(grid),right(grid))
    Mgrid = MaskedGrid(grid,dom)
    plot_grid(Mgrid)
end

function plot_expansion{N,T}(f::FrameFun{N,T}; n=35)
    Tgrid = TensorProductGrid([EquispacedGrid(n, left(set(expansion(f)),idx), right(set(expansion(f)),idx)) for idx = 1:ndims(f)]...)
    data = real(expansion(f)(Tgrid))
    Main.PyPlot.surf(BasisFunctions.range(grid(Tgrid,1)),BasisFunctions.range(grid(Tgrid,2)),data,rstride=1, cstride=1, cmap=Main.PyPlot.ColorMap("coolwarm"),linewidth=0, antialiased=false,vmin=-1.0,vmax=1.0)
end

function plot_boundary(d::AbstractDomain, B::BBox; n =100)
    Tgrid = equispaced_aspect_grid(B,n)
    bound = boundary(Tgrid,d)
    plot_grid(bound)
end
