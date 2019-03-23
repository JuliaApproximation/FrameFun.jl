
"""
Compute a scattered grid of M points randomly distributed in `Ω`, using
the uniform probability measure on `Ω`.
"""
function randomgrid(Ω::Domain, M::Int)
    points = [randompoint(Ω, boundingbox(Ω)) for m in 1:M]
    ScatteredGrid(sort(points))
end

"Generate a single random point inside the given box, with `eltype` `T`."
function randompoint(dom::ProductDomain)
    convert(eltype(dom),map(randompoint,elements(dom)))
end


# Don't return an SVector in 1d, just a value
function randompoint(dom::AbstractInterval)
    T = float(eltype(dom))
    val = rand(T)
    convert(T,val * infimum(dom) + (1-val) * supremum(dom))
end



"""
Generate a single random point inside the given domain, with `eltype` `T`.
Random points are generated inside the given box, until one is inside the domain.
"""
function randompoint(Ω::Domain, box = boundingbox(Ω))
    local point
    indomain = false
    while ~indomain
        point = randompoint(box)
        indomain = point ∈ Ω
    end
    point
end
