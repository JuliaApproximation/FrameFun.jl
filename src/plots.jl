# plots.jl

import BasisFunctions: plot_error

# One-dimensional plot, just the domain
function plot(f::FrameFun{1}; n=201, color = "blue", alpha = 1.0)
    g = grid(resize(basis(f),n))
    mg = MaskedGrid(g, domain(f))
    x = map(Float64, collect(mg))
    data = map(Float64, real(f(mg)))
    Main.PyPlot.plot(x, data, color=color, alpha = alpha)
    Main.PyPlot.title("FrameFun (domain)")
end

function plot_full(f::FrameFun; args...)
    e = SetExpansion(basis(set(f)), coefficients(f))
    BasisFunctions.plot(e; title = "Full extension", args...)
end

function plot_error(f::FrameFun{1}, f_orig::Function; n=201, color="blue")
    g = grid(resize(basis(f),n))
    mg = MaskedGrid(g, domain(f))
    x = collect(mg)
    data = real(f(mg))
    plotdata = abs(Float64[real(f_orig(x)) for x in mg] - data)
    Main.PyPlot.semilogy(x, plotdata, color=color)
    Main.PyPlot.ylim([min(minimum(log10(plotdata)),-16),1])
    Main.PyPlot.title("Absolute Error")
end


function plot_samples(f::FrameFun{1}; gamma=2)
    grid, fbasis2 = oversampled_grid(domain(f), basis(f), gamma)
    x = collect(grid)
    data = real(f(grid))
    Main.PyPlot.stem(x,data)
    Main.PyPlot.title("samples")
end


function plot_domain(d::AbstractDomain2d; n=1000)
    B = boundingbox(d)
    grid = equispaced_aspect_grid(B,n)
    Z = evalgrid(grid, d)
    Main.PyPlot.imshow(Z',interpolation="bicubic",cmap="Blues",extent=(left(B)[1], right(B)[1], left(B)[2], right(B)[2]),aspect="equal",origin="lower")
end



function plot(f::FrameFun{2}; n=101)
    B = boundingbox(domain(set(expansion(f))))
    Tgrid = equispaced_grid(B,n)
    Mgrid = MaskedGrid(Tgrid, domain(f))
    data = map(Float64, real(expansion(f)(Mgrid)))
    pts = collect(Mgrid)
    x = [p[1] for p in pts]
    y = [p[2] for p in pts]
    Main.PyPlot.plot_trisurf(x, y, data)
end

function plot_image(f::FrameFun{2}; n=301,unscaled=false, border=true, colorbar=true)
    d = domain(f)
    Tgrid = grid(resize(basis(f), (n,n)))
    Mgrid = MaskedGrid(Tgrid, d)
    data = real(f(Tgrid))
    pts = collect(Tgrid)
    x = [p[1] for p in pts]
    y = [p[2] for p in pts]
    X = reshape(x, size(Tgrid))
    Y = reshape(y, size(Tgrid))

    for i in eachindex(Tgrid)
        if ~in(i, Mgrid)
            data[i] = NaN
        end
    end

    my_cmap = Main.PyPlot.matplotlib[:cm][:get_cmap]("Spectral",100)
    my_cmap[:set_bad](color="#663300", alpha=0)
    my_cmap[:set_under](color="#663300", alpha=0)
    Main.PyPlot.plt[:register_cmap](name="my_cmap",cmap=my_cmap)

    vmin = minimum(data)
    vmax = maximum(data)
    Main.PyPlot.pcolormesh(X, Y, data, vmin=vmin, vmax=vmax, shading = "gouraud",
        cmap = my_cmap)
    # bound = boundary(Tgrid,d)
    # border && plot_grid(bound)
    Main.PyPlot.axis("scaled")
    Main.PyPlot.xlim([left(Tgrid)[1], right(Tgrid)[1]])
    Main.PyPlot.ylim([left(Tgrid)[2], right(Tgrid)[2]])
    colorbar && Main.PyPlot.colorbar()
end

function plot_error(f::FrameFun{2}, f_orig::Function; n=301, border=false)
    d = domain(f)
    Tgrid = grid(resize(basis(f), (n,n)))
    Mgrid = MaskedGrid(Tgrid, d)
    f_data = real(f(Tgrid))
    orig_data = Float64[f_orig(x...) for x in Tgrid]
    orig_data = reshape(orig_data, size(Tgrid))
    data = f_data - orig_data;
    pts = collect(Tgrid)
    x = [p[1] for p in pts]
    y = [p[2] for p in pts]
    X = reshape(x, size(Tgrid))
    Y = reshape(y, size(Tgrid))

    for i in eachindex(Tgrid)
        if ~in(i, Mgrid)
            data[i] = NaN
        end
    end

    my_cmap = Main.PyPlot.matplotlib[:cm][:get_cmap]("Spectral",100)
    my_cmap[:set_bad](color="#663300", alpha=0)
    my_cmap[:set_under](color="#663300", alpha=0)
    Main.PyPlot.plt[:register_cmap](name="my_cmap",cmap=my_cmap)

    vmin = minimum(data)
    vmax = maximum(data)
    Main.PyPlot.pcolormesh(X,Y,data,vmin=vmin, vmax=vmax, shading = "gouraud",
        cmap = my_cmap)
    # bound = boundary(Tgrid,d)
    # border && plot_grid(bound)
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
