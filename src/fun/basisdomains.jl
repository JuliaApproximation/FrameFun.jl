# basisdomains.jl

# We associate a domain with some of the bases defined in BasisFunctions.

domain(b::FourierBasis) = interval(0,1)

domain(b::ChebyshevBasis) = interval(-1, 1)

domain(d::TensorProductDict) = cartesianproduct(map(domain, elements(d))...)

domain(dict::MappedDict) = apply_map(domain(superdict(dict)), mapping(dict))
