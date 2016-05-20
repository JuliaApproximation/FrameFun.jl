# explicit_extension.jl

immutable LinearExtension{G <: AbstractSubGrid,ELT} <: AbstractOperator{ELT}
    grid    ::  G
end

LinearExtension{G <: AbstractSubGrid,T}(s::DiscreteGridSpace{G,1,T}) = LinearExtension{G,T}(grid(s))

src(op::LinearExtension) = DiscreteGridSpace(op.grid, eltype(op))

dest(op::LinearExtension) = DiscreteGridSpace(grid(op.grid), eltype(op))


apply!(op::LinearExtension, coef_dest, coef_src) =
    apply_extension!(op, op.grid, grid(op.grid), coef_dest, coef_src)

function apply_extension!(op::LinearExtension, subgrid::IndexSubGrid, supergrid::AbstractGrid, coef_dest, coef_src)
    i1 = subgrid.i1
    i2 = subgrid.i2
    t1 = supergrid[i1]
    t2 = supergrid[i2]
    val1 = coef_src[1]
    val2 = coef_src[end]
    dist = supergrid[i1] - left(supergrid) + right(supergrid) - supergrid[i2]
    for i in i1:i2
        coef_dest[i] = coef_src[i-i1+1]
    end
    for i in 1:i1-1
        coef_dest[i] = val1 + (t1 - supergrid[i])/dist * (val2-val1)
    end
    for i in i2+1:length(supergrid)
        coef_dest[i] = val2 - (supergrid[i]-t2)/dist * (val2-val1)
    end
    coef_dest
end
