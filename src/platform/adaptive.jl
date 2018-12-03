
abstract type CriteriumStyle end

struct ResidualStyle <: CriteriumStyle end


# Generic adaptivity
function approximate(platform::Platform, fun; options...)
end
