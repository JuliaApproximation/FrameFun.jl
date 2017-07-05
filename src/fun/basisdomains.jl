# basisdomains.jl

# We associate a domain with some of the bases defined in BasisFunctions.

domain(b::FourierBasis) = interval(0,1)

domain(b::ChebyshevBasis) = interval(-1, 1)

domain(s::TensorProductSet) = tensorproduct(map(domain, elements(s))...)

domain(set::MappedSet) = apply_map(domain(superset(set)), mapping(set))
