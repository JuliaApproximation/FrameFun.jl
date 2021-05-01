
export Characteristic
"""
A domain that is described by its characteristic function.
"""
struct Characteristic{N,T} <: Domain{SVector{N,T}}
    char    ::  Function
    box     ::  Domain{SVector{N,T}}
end

export characteristic
characteristic(char::Function, dom)  = Characteristic(char,boundingbox(dom))

DomainSets.indomain(x, c::Characteristic) = c.char(x)

boundingbox(c::Characteristic) = c.box

show(io::IO, c::Characteristic) = print(io, "a domain described by a characteristic function")
