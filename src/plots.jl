# plots.jl

# One-dimensional plot, just the domain
function plot(f::SetExpansion;n=200)
    grid=EquispacedGrid(n,left(domain(f)),right(domain(f)))
    data=f(grid)
    Main.PyPlot.plot(BasisFunctions.range(grid),data)
end

## # One-dimensional plot, including extension
## function plot_full{1}(f::Fun{1})
    
## end
function plot_expansion(f::SetExpansion;n=200)
    grid=EquispacedGrid(n,left(basis(set(f))),right(basis(set(f))))
    data = f(grid)
    Main.PyPlot.plot(BasisFunctions.range(grid),data)
end

# Maybe place this in funs.jl?
function call(f::FrameFun, g::AbstractGrid)
    result = Array(eltype(f), size(g))
    call!(f, result, g)
    result
end

function call!{N}(f::FrameFun, result::AbstractArray, g::AbstractGrid{N})
    x = Array(eltype(f), N)
    for i in eachindex(g)
        getindex!(x, g, i)
        result[i] = call(f, x...)
    end
end

function call!(f::FrameFun, result::AbstractArray, g::AbstractGrid1d)
    for i in eachindex(g)
        result[i] = call(f, g[i])
    end
end

function call!(f::FrameFun, result::AbstractArray, x::AbstractArray)
    @assert size(result) == size(x)
    for i = 1:length(x)
        result[i] = call(f, x[i])
    end
end



## function plot{N,T}(f::FrameFun{N,T};n=100)
##     Tgrid=TensorProductGrid([EquispacedGrid(100, left(set(expansion(f)),idx), right(set(expansion(f)),idx)) for idx = 1:N]...)
##     Mgrid=MaskedGrid(Tgrid, domain(f))
##     data = expansion(f)(Mgrid)
##     x=[Mgrid[i][1] for i = 1:length(Mgrid)]
##     y=[Mgrid[i][2] for i = 1:length(Mgrid)]
##     Main.PyPlot.plot_trisurf(x,y,data)
## end

function plot{N,T}(f::FrameFun{N,T};n=35)
    Tgrid=TensorProductGrid([EquispacedGrid(n, left(box(domain(f)),idx), right(box(domain(f)),idx)) for idx = 1:N]...)
    data = real(f(Tgrid))
    Main.PyPlot.surf(BasisFunctions.range(grid(Tgrid,1)),BasisFunctions.range(grid(Tgrid,2)),data,rstride=1, cstride=1, cmap=Main.PyPlot.ColorMap("coolwarm"),linewidth=0, antialiased=false,vmin=-1.0,vmax=1.0)
end
## # One-dimensional plot, including extension
## function plot_full{1}(f::Fun{1})
    
## end
function plot_expansion{N,T}(f::FrameFun{N,T};n=35)
    Tgrid=TensorProductGrid([EquispacedGrid(n, left(set(expansion(f)),idx), right(set(expansion(f)),idx)) for idx = 1:N]...)
    data = real(expansion(f)(Tgrid))
    Main.PyPlot.surf(BasisFunctions.range(grid(Tgrid,1)),BasisFunctions.range(grid(Tgrid,2)),data,rstride=1, cstride=1, cmap=Main.PyPlot.ColorMap("coolwarm"),linewidth=0, antialiased=false,vmin=-1.0,vmax=1.0)
end

