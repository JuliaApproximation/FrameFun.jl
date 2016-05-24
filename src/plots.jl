# plots.jl

import BasisFunctions: plot_error

# One-dimensional plot, just the domain
function plot(f::FrameFun1d; n=201, color = "blue", alpha = 1.0)
    vals, x, ~ = evaluate_framefun(f, n)

    x = map(Float64, x)
    vals = map(Float64, vals)

    plot(x, vals, color=color, alpha = alpha)
    xlim([left(domain(f)), right(domain(f))])
    title("FrameFun (domain)")
end

function evaluate_framefun(f::FrameFun1d, n; full = false)
    fgrid = grid(resize(basis(f),n))
    mgrid = subgrid(fgrid, domain(f))
    x = collect(fgrid)
    vals = real(f(fgrid))

    if ~full
        for i in eachindex(fgrid)
            if ~in(i, mgrid)
                vals[i] = NaN
            end
        end
    end
    vals, x, fgrid
end

function evaluate_framefun(f::FrameFun2d, n; full = false)
    fgrid = grid(resize(basis(f), (n,n)))
    mgrid = subgrid(fgrid, domain(f))
    vals = real(f(fgrid))
    pts = collect(fgrid)
    x = [p[1] for p in pts]
    y = [p[2] for p in pts]
    X = reshape(x, size(fgrid))
    Y = reshape(y, size(fgrid))

    if ~full
        for i in eachindex(fgrid)
            if ~in(i, mgrid)
                vals[i] = NaN
            end
        end
    end
    vals, X, Y, fgrid
end


function plot_extension(f::FrameFun; args...)
    e = SetExpansion(basis(set(f)), coefficients(f))
    plot(e; title = "Full extension", args...)
end

function plot_error(f::FrameFun1d, f_orig::Function; n=201, color="blue")
    vals, x, fgrid = evaluate_framefun(f, n)
    errorvals = abs(eltype(f)[real(f_orig(x)) for x in fgrid] - vals)

    x = map(Float64, x)
    errorvals = map(Float64, errorvals)

    semilogy(x, errorvals, color=color)
    xlim([left(domain(f)), right(domain(f))])
    ylim([min(minimum(log10(errorvals)),-16),1])
    title("Absolute Error")
end


function plot_samples(f::FrameFun{1}; gamma=2)
    samplegrid,~ = oversampled_grid(domain(f), basis(f), gamma)
    x = collect(samplegrid)
    vals = real(f(samplegrid))
    stem(x, vals)
    title("samples")
end


function plot(d::AbstractDomain2d; n=1000)
    B = boundingbox(d)
    grid = equispaced_aspect_grid(B,n)
    Z = evalgrid(grid, d)
    imshow(Z',interpolation="bicubic",cmap="Blues",extent=(left(B)[1], right(B)[1], left(B)[2], right(B)[2]),aspect="equal",origin="lower")
end


function plot_surf(f::FrameFun{2}; n=101)
    B = boundingbox(domain(set(expansion(f))))
    fgrid = equispaced_grid(B,n)
    Mgrid = MaskedGrid(fgrid, domain(f))
    data = map(Float64, real(expansion(f)(Mgrid)))
    pts = collect(Mgrid)
    x = [p[1] for p in pts]
    y = [p[2] for p in pts]
    plot_trisurf(x, y, data)
end

function plot(f::FrameFun{2}; n=301, colorbar=true)
    vals, X, Y, fgrid = evaluate_framefun(f, n)

    my_cmap = matplotlib[:cm][:get_cmap]("Spectral",100)
    my_cmap[:set_bad](color="#663300", alpha=0)
    my_cmap[:set_under](color="#663300", alpha=0)
    plt[:register_cmap](name="my_cmap",cmap=my_cmap)

    vmin = minimum(vals)
    vmax = maximum(vals)
    pcolormesh(X, Y, vals, vmin=vmin, vmax=vmax, shading = "gouraud",
        cmap = my_cmap)
    axis("scaled")
    xlim([left(fgrid)[1], right(fgrid)[1]])
    ylim([left(fgrid)[2], right(fgrid)[2]])
    colorbar && PyPlot.colorbar()
end

function plot_error(f::FrameFun{2}, f_orig::Function; n=301, colorbar = true)
    vals, X, Y, fgrid = evaluate_framefun(f, n, full = true)
    orig_vals = Float64[f_orig(x...) for x in fgrid]
    orig_vals = reshape(orig_vals, size(fgrid))
    data = log10(abs(vals - orig_vals))

    # my_cmap = matplotlib[:cm][:get_cmap]("Spectral",100)
    # my_cmap[:set_bad](color="#663300", alpha=0)
    # my_cmap[:set_under](color="#663300", alpha=0)
    # plt[:register_cmap](name="my_cmap",cmap=my_cmap)

    vmin = minimum(data)
    vmax = maximum(data)
    pcolormesh(X, Y, data, vmin=vmin, vmax=vmax, shading = "gouraud")
    axis("scaled")
    xlim([left(fgrid)[1], right(fgrid)[1]])
    ylim([left(fgrid)[2], right(fgrid)[2]])
    bound = boundary(fgrid, domain(f))
    plot_grid(bound)
    colorbar && PyPlot.colorbar()
end

function plot_grid(grid::AbstractGrid2d)
    pts = collect(grid)
    x = [p[1] for p in pts]
    y = [p[2] for p in pts]
    plot(x,y,linestyle="none",marker="o",color="black",mec="black",markersize=1.5)
    axis("equal")
end

function plot_grid(grid::AbstractGrid3d)
    pts = collect(grid)
    x = [p[1] for p in pts]
    y = [p[2] for p in pts]
    z = [p[3] for p in pts]
    plot3D(x,y,z,linestyle="none",marker="o",color="blue")
    axis("equal")
end

function plot_boundary(d::AbstractDomain, B::BBox; n =100)
    fgrid = equispaced_aspect_grid(B,n)
    bound = boundary(fgrid,d)
    plot_grid(bound)
end
