# A collection of extensions to the DomainSets package.

using DomainSets: inverse_map, forward_map


##########################################################################
### Distances and Normals
##########################################################################
export normal, distance, volume
# Function that returns the normal vector when evaluated at the domain boundary
# This vector is of unit length and points to the OUTSIDE.
normal(x, d::Domain) = error("Normal not available for this domain type")

# Auxiliary function that returns distance from the boundary in some metric
distance(x, d::Domain) = error("Domain distance not available for this domain type")

distance(x, ::UnitSimplex) = min(minimum(x),1-sum(x))

function normal(x, ::UnitSimplex)
    z = fill(eltype(x)(0), size(x))
    if minimum(x)<abs(sum(x)-1)/sqrt(length(x))
        index = findmin(x)[2]
        setindex!(z,-1,index)
    else
        z = fill(eltype(x)(1), size(x))
    end
    return z/norm(z)
end

normal(x, d::UnitHyperBall) = x/norm(x)

distance(x, d::AbstractInterval) = min(maximum(d)-x,x-minimum(d))

normal(x, d::AbstractInterval) = (abs(minimum(d)-x) < abs(maximum(d)-x)) ? -1 : 1

distance(x, d::UnitHyperBall) = 1-norm(x)

normal(x, d::UnitHyperSphere) = x/norm(x)

distance(x, d::UnitHyperSphere) = 1-norm(x)

distance(x,d::UnionDomain) = indomain(x,d) ? sum(map(di->max(0,distance(x,di)),elements(d))) : maximum(map(di->distance(x,di),elements(d)))

normal(x,d::UnionDomain) = normal(x,elements(d)[findmin(map(di->abs(distance(x,di)),elements(d)))[2]])

distance(x,d::IntersectionDomain) = minimum(map(di->distance(x,di),elements(d)))

normal(x,d::IntersectionDomain) = normal(x,elements(d)[findmin(map(di->distance(x,di),elements(d)))[2]])

distance(x,d::DifferenceDomain) = indomain(x,d) ? min(abs(distance(x,d.domains[1])),abs(distance(x,d.domains[2]))) : -1*min(abs(distance(x,d.domains[1])),abs(distance(x,d.domains[2])))

normal(x,d::DifferenceDomain) = abs(distance(x,d.domains[1]))<abs(distance(x,d.domains[2])) ? normal(x,d.domains[1]) : -1*normal(x,d.domains[2])

distance(x, t::ProductDomain) = minimum(map(distance, x, elements(t)))

distance(x, d::DomainSets.MappedDomain) = distance(inverse_map(d)*x,source(d))

function normal(x, d::DomainSets.MappedDomain)
    x = applymap(inverse_map(d),normal(inverse_map(d)*x,source(d)))
    x0 = apply_inverse(inverse_map(d),zeros(size(x)))
   (x-x0)/norm(x-x0)
end
function normal(x, t::ProductDomain)
    index = findmin(map(distance, x, elements(t)))[2]
    [(i==index)*normal(x[i],element(t,i)) for i =1:length(elements(t))]
end


##########################################################################
### Volume
##########################################################################

volume(d::Domain) = missing

volume(d::AbstractInterval) = abs(rightendpoint(d)-leftendpoint(d))


# TODO: remove this again in the future
DomainSets.approx_indomain(x, d::DifferenceDomain, tolerance) = in(x, d)
