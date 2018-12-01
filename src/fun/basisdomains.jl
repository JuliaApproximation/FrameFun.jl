
# We associate a domain with some of the bases defined in BasisFunctions.

# TODO: remove now that BasisFunctions has `support` implemented
domain(d::TensorProductDict) = cartesianproduct(map(domain, elements(d))...)

domain(dict::MappedDict) = apply_map(domain(superdict(dict)), mapping(dict))
