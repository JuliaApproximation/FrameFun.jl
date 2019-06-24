export param_first, param_next, param_increment, param_inbetween
"""
    param_first(platform::Platform)
"""
param_first(platform::Platform) = 1
"""
    param_next(platform::Platform, param)
"""
param_next(platform::Platform, param) =
    extension_size(dictionary(platform, param))
"""
    param_increment(platform::Platform, param)
"""
param_increment(platform::Platform, n) = addone(n)
"""
    param_inbetween(platform::Platform, param1, param2)
"""
param_inbetween(platform::Platform, param1, param2) =
    inbetween(param1, param2)

addone(x::Number) = x+1
addone(x) = map(addone, x)
inbetween(n1::Int, n2::Int) = (n1+n2) >> 1
inbetween(n1::T, n2::T) where {T} = map(inbetween, n1, n2)
