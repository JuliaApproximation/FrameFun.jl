
################################################################################
### A polar domain described by a single variable function
################################################################################
export PolarDomain
struct PolarDomain{T} <: Domain{SVector{2,T}}
    charFun
    box     ::  Domain{SVector{2,T}}
end
export polardomain
function polardomain(char, domain)
    T = subeltype(domain)
    # F = Fun(char, Fourier(100, -T(pi), T(pi)), Interval(-T(pi), T(pi)))
    # PolarDomain(F, boundingbox(domain))
    PolarDomain(char, boundingbox(domain))
end

indomain(x, c::PolarDomain) = sqrt(x[1]^2+x[2]^2) < real(c.charFun(atan(x[2],x[1])))

boundingbox(c::PolarDomain) = c.box

show(io::IO, c::PolarDomain) = print(io, "a domain in polar coordinates")

distance(x, c::PolarDomain) = real(c.charFun(atan(x[2],x[1])))-norm(x)

function normal(c::PolarDomain, x)
    phi = atan(x[2], x[1])
    y=c.charFun'(phi)*sin(phi)+c.charFun(phi)*cos(phi)
    x=c.charFun'(phi)*cos(phi)-c.charFun(phi)*sin(phi)
    return [real(y/(norm([x;y]))),-real(x/(norm([x;y])))]
end
