
"""
Makes sure that (i-1)/N < x < i/N holds
If x ≈ i/N it return -i
"""
interval_index(B::BSplineTranslatesBasis,x::Real) = round(x*length(B))≈x*length(B) ? -round(Int,x*length(B))-1 : ceil(Int,x*length(B))

function BasisFunctions.support(B::BSplineTranslatesBasis, i)
    start = (i-1)/length(B)
    width = (degree(B)+1)/length(B)
    stop  = start+width
    stop <=1 ? (return [start,stop]) : (return [0.,stop-1.], [start,1.])
end

# convert 2D indices to linear index
# TODO is it possible to generalize this function to ND
function create_indices(B, i1, i2)
    [linear_index(B,(i,j)) for i in i1 for j in i2]
#     [(i,j) for i in i1 for j in i2]
end

function create_indices(B, i1, i2, i3)
    [linear_index(B,(i,j,k)) for i in i1 for j in i2 for k in i3]
end

"""
The linear index of the spline elements of B that are non-zero in x.
"""
function overlapping_spline_indices(B::BSplineTranslatesBasis, x::Real)
    # The interval_index is the starting index of all spline elements that overlap with x
    init_index = interval_index(B,x)
    degree(B) == 0 && return abs(init_index)
    # The number of elements that overlap with one interval
    no_elements = degree(B)+1
    if init_index < 0
        init_index = -init_index-1
        no_elements = no_elements-1
    end
    [mod(init_index+i-2,length(B)) + 1 for i in 1:-1:2-no_elements]
end

function overlapping_spline_indices(B::TensorProductDict, x::SVector)
    index_sets = [overlapping_spline_indices(s,x[i]) for (i,s) in enumerate(elements(B))]
    create_indices(B,index_sets...)
end

"""
Index of elements of `B` that overlap with `boundary`.
"""
function boundary_element_indices(B, boundary::MaskedGrid)
    s = Set{Int}()
    for x in boundary
        push!(s,overlapping_spline_indices(B,x)...)
    end
    collect(s)
end

"""
The linear indices of the points of `g` at which B[i] is not zero.
"""
function support_indices(B::BSplineTranslatesBasis, g::AbstractEquispacedGrid, i)
    indices = Vector{Int}()
    dx = stepsize(g)
    x0 = g[1]
    s = support(B,i)
    if length(s[1]) == 1
        start = ceil(Int,(s[1]-x0)/dx)
        stop = floor(Int,(s[2]-x0)/dx)
        (degree(B) != 0) && ((s[1]-x0)/dx ≈ start) && (start += 1)
        ((s[2]-x0)/dx ≈ stop) && (stop -= 1)
        push!(indices,(start+1:stop+1)...)
    else
        interval = s[1]
        start = 0
        stop = floor(Int,(interval[2]-x0)/dx)
        ((interval[2]-x0)/dx ≈ stop) && (stop -= 1)
        push!(indices,(start+1:stop+1)...)

        interval = s[2]
        start = ceil(Int,(interval[1]-x0)/dx)
        stop = length(g)-1
        ((interval[1]-x0)/dx ≈ start) && (start += 1)
        push!(indices,(start+1:stop+1)...)
    end
end

function support_indices(B::TensorProductDict, g::ProductGrid, index::Int)
    cartindex = ind2sub(size(B),index)
    index_sets = [support_indices(s,element(g,i),cartindex[i]) for (i,s) in enumerate(elements(B))]
    create_indices(g,index_sets...)
end

"""
A grid that contains the points of `omega_grid` that are not evaluated to zero by the elements that overlap with boundary_grid.
"""
function boundary_support_grid(B, boundary_grid::MaskedGrid, omega_grid::MaskedGrid)
    boundary_indices = boundary_element_indices(B,boundary_grid)
    s = Set{Int}()
    for i in boundary_indices
        push!(s,support_indices(B,supergrid(omega_grid),i)...)
    end
    a = collect(s)
    mask = zeros(Bool,size(omega_grid.mask))
    mask[a] = true
    mask .= mask .& omega_grid.mask
    MaskedGrid(supergrid(omega_grid),mask)
end

"""
An index extension operator from the elements of `B` the overlap with the boundary to the complete dictionary.
"""
function boundary_extension_operator(boundary::AbstractGrid, B::Dictionary)
    boundary_indices = boundary_element_indices(B,boundary)
    IndexExtensionOperator(Span(B[boundary_indices]),Span(B), boundary_indices)
end
