# randomgrid.jl

"""
Compute a scattered grid of M points randomly distributed in `Ω`, using
the uniform probability measure on `Ω`.
"""
randomgrid(Ω::Domain, M::Int) =
    ScatteredGrid([randompoint(Ω, boundingbox(Ω)) for m in 1:M])

"Generate a single random point inside the given box, with `eltype` `T`."
function randompoint(dom::ProductDomain) 
    convert(eltype(dom),map(randompoint,elements(dom)))
end

# Don't return an SVector in 1d, just a value
function randompoint(dom::Interval)
    # Random is not implemented for BigFloats, see issue #13950
    # That's why we just use Float64
    val = rand()
    convert(eltype(dom),val * leftendpoint(dom) + (1-val) * rightendpoint(dom))
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
