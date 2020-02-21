export param_first, param_double, param_increment, param_inbetween, hasparam_inbetween
"""
    param_first(platform::Platform)
"""
param_first(platform::Platform) = 1
"""
    param_double(platform::Platform, param)
"""
param_double(platform::Platform, param) =
    extensionsize(dictionary(platform, param))
"""
    param_increment(platform::Platform, param)
"""
param_increment(platform::Platform, n) = addone(n)
"""
    param_inbetween(platform::Platform, param1, param2)
"""
param_inbetween(platform::Platform, param1, param2) =
    inbetween(param1, param2)
"""
    hasparam_inbetween(platform::Platform, param1, param2, stoptolerance)
"""
hasparam_inbetween(platform::Platform, param1, param2, stoptolerance::Real) =
    hasparam_inbetween(param1, param2, stoptolerance)

addone(x::Number) = x+1
addone(x) = map(addone, x)
inbetween(n1::Int, n2::Int) = (n1+n2) >> 1
hasparam_inbetween(n1::Int, n2::Int, stoptolerance) = inbetween(n1, n2) != n1
inbetween(n1::T, n2::T) where {T<:Tuple} = map(inbetween, n1, n2)
hasparam_inbetween(n1::T, n2::T, stoptolerance) where {T<:Tuple} =
    any(map(hasparam_inbetween, n1, n2, Ref(stoptolerance))...)
