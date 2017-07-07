# characteristic.jl

"""
A domain that is described by its characteristic function.
"""
struct Characteristic{N,T} <: EuclideanDomain{N,T}
    char    ::  Function
    box     ::  EuclideanDomain{N,T}
end

characteristic(char::Function, dom)  = Characteristic(char,boundingbox(dom))

indomain(x, c::Characteristic) = c.char(x)

boundingbox(c::Characteristic) = c.box

show(io::IO, c::Characteristic) = print(io, "a domain described by a characteristic function")
