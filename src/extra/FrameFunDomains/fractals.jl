
###############################################################################
## The Mandelbrot set
###############################################################################
export Mandelbrot,mandelbrot
# TODO: fractals should be subsets of the complex plane
struct Mandelbrot{T} <: Domain{SVector{2,T}}
    maxiter     ::  Int
    threshold   ::  T
    maskcache   ::  Dict
    box         ::  Domain
end

Base.isempty(::Mandelbrot) = false

function mandelbrot(maxiter = 1000, threshold = 1000.0)
    box = rectangle(-1.0, 0.35, -0.65, 0.65)
    M1 = 136
    M2 = 200
    mask1 = computemandelbrotgrid(EquispacedGrid(M1,-1.0, 0.3)× EquispacedGrid(M2,-0.65, 0.65), maxiter, threshold)
    mask2 = computemandelbrotgrid(EquispacedGrid(M1,-1.0, 0.3)× EquispacedGrid(M2,-0.65, 0.65), maxiter, threshold)
    cache = Dict{Int,BitArray{2}}()
    cache[M1] = mask1
    cache[M2] = mask2
    Mandelbrot(maxiter, threshold, cache, box)
end


function mandelbrotiteration(x, maxiter, threshold)
    c = 2*(x[1]+1im*x[2])
    T = typeof(c)
    z = zero(T)
    iter = 0
    while abs(z) < threshold && iter < maxiter
        z = z^2 + c
        iter += 1
    end
    abs(z) < threshold
end

function computemandelbrotgrid(grid, maxiter, threshold)
    mask = BitArray(undef, size(grid))
    fill!(mask,0)
    for i_2 = 1:size(grid, 2)
        for i_1 = 1:size(grid, 1)
            mask[i_1,i_2] = mandelbrotiteration(grid[i_1,i_2], maxiter, threshold)
        end
    end
    mask
end

# TODO: reconsider broadcast
function indomain_broadcast(grid::ProductGrid, m::Mandelbrot)
    if (infimum(support(grid)) ≈ minimum(boundingbox(m))) && (supremum(support(grid)) ≈ maximum(boundingbox(m)))
        if haskey(m.maskcache, size(grid,1))
            mask = m.maskcache[size(grid,1)]
        else # compute mask and cache it
            mask = computemandelbrotgrid(grid, m.maxiter, m.threshold)
            m.maskcache[size(grid,1)] = mask
        end
    else # Don't cache if the grid doesn't match the bounding box
        mask = computemandelbrotgrid(grid, m.maxiter, m.threshold)
    end
    mask
end

function isapprox(t::Tuple{T,T}, v::SVector{2,T}) where {T}
    return t[1]≈v[1] && t[2]≈v[2]
end
function isapprox(v::SVector{2,T}, t::Tuple{T,T}) where {T}
    return t[1]≈v[1] && t[2]≈v[2]
end
boundingbox(m::Mandelbrot) = m.box

indomain(x::SVector{2}, m::Mandelbrot) = mandelbrotiteration(x, m.maxiter, m.threshold)

show(io::IO, m::Mandelbrot) = print(io, "The Mandelbrot set")


################
## Julia sets
################
export JuliaSet, juliaset
struct JuliaSet{T} <: Domain{SVector{2,T}}
    c           ::  Complex{T}
    maxiter     ::  Int
    maskcache   ::  Dict
    box         ::  Domain
end

Base.isempty(::JuliaSet) = false

function juliaset(c = -0.122565+0.744866im, maxiter = 1000)
    box = rectangle(-0.2, 1.2, -0.4, 0.4)

    mask1 = computejuliasetgrid(EquispacedGrid(200,-0.2, 1.2)× EquispacedGrid(200,-0.4, 0.4), c, maxiter)
    mask2 = computejuliasetgrid(EquispacedGrid(200,-0.2, 1.2)× EquispacedGrid(200,-0.4, 0.4), c, maxiter)
    cache = Dict{Int,BitArray{2}}()
    cache[100] = mask1
    cache[200] = mask2
    JuliaSet(c, maxiter, cache, box)
end


function juliasetiteration(x, c, maxiter)
    gamma = 1 + sqrt(1-4*c)
    z = x[1] + 1im*x[2]
    for i = 1:maxiter
        z = gamma*z*(1-z)
    end
    abs(z) < 1000
end

function computejuliasetgrid(grid, c, maxiter)
    m = size(grid)
    mask = BitArray(undef, m)
    fill!(mask, 0)
    for i_2 = 1:m[2]
        for i_1 = 1:m[1]
            a,b = grid[i_1,i_2]
            mask[i_1,i_2] = juliasetiteration(grid[i_1,i_2], c, maxiter)
        end
    end
    mask
end

# TODO: reconsider broacast
function indomain_broadcast(grid, js::JuliaSet)
    if (map(x->infimum(support(x)), components(grid)) ≈ infimum(boundingbox(js))) && (map(x->supremum(support(x)),components(grid)) ≈ supremum(boundingbox(js)))
        if haskey(js.maskcache, size(grid,1))
            mask = js.maskcache[size(grid,1)]
        else # compute mask and cache it
            mask = computejuliasetgrid(grid, js.c, js.maxiter)
            js.maskcache[size(grid,1)] = mask
        end
    else # Don't cache if the grid doesn't match the bounding box
        mask = computejuliasetgrid(grid, js.c, js.maxiter)
    end
    mask
end

indomain(x::SVector{2}, js::JuliaSet) = juliasetiteration(x, js.c, js.maxiter)

show(io::IO, js::JuliaSet) = print(io, "A particular Julia Set also known as the Douady rabbit")

boundingbox(js::JuliaSet) = js.box
