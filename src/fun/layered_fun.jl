# layered_fun.jl

"""
A LayeredFun consists of several layers of functions, that can be on top of each other.
"""
structLayeredFun
    domains     ::  Vector{AbstractDomain}
    funs        ::  Vector{AbstractFun}
end
