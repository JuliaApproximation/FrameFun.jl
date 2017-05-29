# curve.jl

abstract type ParametricCurve{N} end

export Circle, point

struct Circle{T} <: ParametricCurve{2}
    radius  ::  T
end

Circle() = Circle(1)

parameter_domain(c::Circle) = Interval(0, 2*π)

point(c::Circle, t) = SVector{2}(c.radius*cos(t), c.radius*sin(t))
