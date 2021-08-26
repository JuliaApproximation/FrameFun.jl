# A collection of extensions to the DomainSets package.

using DomainSets: inverse_map, forward_map, MappedDomain

##########################################################################
### Distances and Normals
##########################################################################
export normal, distance
# Function that returns the normal vector when evaluated at the domain boundary
# This vector is of unit length and points to the OUTSIDE.
function normal(d::Domain, x)
    if hascanonicaldomain(d)
        normal(convert(MappedDomain, d), x)
    else
        error("Normal not available for this domain type")
    end
end

# Auxiliary function that returns distance from the boundary in some metric
function distance(x, d::Domain)
    if hascanonicaldomain(d)
        distance(x, convert(MappedDomain, d))
    else
        error("Domain distance not available for this domain type")
    end
end

distance(x, ::UnitSimplex) = min(minimum(x),1-sum(x))

distance(x, d::AbstractInterval) = min(abs(maximum(d)-x),abs(x-minimum(d)))

distance(x, d::UnitBall) = 1-norm(x)

distance(x, d::UnitSphere) = 1-norm(x)

distance(x,d::UnionDomain) = indomain(x,d) ? sum(map(di->max(0,distance(x,di)),components(d))) : maximum(map(di->distance(x,di),components(d)))

normal(d::UnionDomain, x) = normal(components(d)[findmin(map(di->abs(distance(x,di)),components(d)))[2]], x)

distance(x,d::IntersectDomain) = minimum(map(di->distance(x,di),components(d)))

normal(d::IntersectDomain, x) = normal(components(d)[findmin(map(di->distance(x,di),components(d)))[2]], x)

distance(x,d::SetdiffDomain) = indomain(x,d) ? min(abs(distance(x,d.domains[1])),abs(distance(x,d.domains[2]))) : -1*min(abs(distance(x,d.domains[1])),abs(distance(x,d.domains[2])))

normal(d::SetdiffDomain, x) = abs(distance(x,d.domains[1]))<abs(distance(x,d.domains[2])) ? normal(d.domains[1], x) : -1*normal(d.domains[2], x)

distance(x, t::ProductDomain) = minimum(map(distance, x, components(t)))

distance(x, d::MappedDomain) = distance(inverse_map(d)(x),superdomain(d))

function normal(d::MappedDomain, x)
    x = applymap(forward_map(d),normal(superdomain(d),inverse_map(d)(x)))
    x0 = inverse(inverse_map(d))(zeros(size(x)))
   (x-x0)/norm(x-x0)
end
function normal(t::ProductDomain, x)
    index = findmin(map(distance, x, components(t)))[2]
    [(i==index)*normal(component(t,i),x[i]) for i =1:length(components(t))]
end
