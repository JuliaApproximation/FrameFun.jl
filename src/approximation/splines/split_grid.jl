using AbstractTrees
import AbstractTrees: children, printnode
import Base: start, done, next
import BasisFunctions: grid

"The node of a tree containing the domain composition of a grid"
abstract type AbstractDomainDecompositionNode end

"""
The leaf of a domain decomposition tree
"""
struct DomainDecompositionLeaf <: AbstractDomainDecompositionNode
    "The grid in which the error is minimized"
    DMZ::Union{BasisFunctions.AbstractGrid}
    "The elemens overlapping this grid can be used to minimize the error"
    grid::Union{Void,BasisFunctions.AbstractGrid}
end

BasisFunctions.grid(leaf::DomainDecompositionLeaf) = leaf.grid
DMZ(leaf::DomainDecompositionLeaf) = leaf.DMZ

"""
A container for multiple decomposition decomposition nodes
"""
mutable struct DomainDecompositionNodeContainer
    container::Vector{N} where {N<: AbstractDomainDecompositionNode}
end

DomainDecompositionNodeContainer(DMZs::Vector{G}) where {G<:AbstractGrid} =
    DomainDecompositionNodeContainer(map(dmz->DomainDecompositionLeaf(dmz, nothing), DMZs))
DomainDecompositionNodeContainer(DMZs::Vector{G1}, grids::Vector{G2}) where {G1<:AbstractGrid, G2<:AbstractGrid} =
    DomainDecompositionNodeContainer([DomainDecompositionLeaf(DMZs[i], grids[i]) for i in 1:length(grids)])
DomainDecompositionNodeContainer(DMZ::AbstractGrid, grid::Union{Void,AbstractGrid}) = DomainDecompositionNodeContainer([DomainDecompositionLeaf(DMZ, grid)])

"""
A node of a domain decomposition tree.

splits in two parts:
    Mid
    Other
Mid is a group of grids which separate the grids in Other.

More details on the separation:
The width of the grids in Mid is such that there is a space the size of the span
of a basis element in which the coefficients can be solved exactly as though they
were solved in the grid that is the combination of Mid and Other.
Changing the problem in on of the Other nodes does not effect the problem each
of the other Other nodes.
"""
struct DomainDecompositionNode <: AbstractDomainDecompositionNode
    Mid ::DomainDecompositionNodeContainer
    Other::DomainDecompositionNodeContainer
end

"Depth of the tree. A leaf is depth 0."
depth(::DomainDecompositionLeaf) = 0
depth(c::DomainDecompositionNodeContainer) = length(c.container)==0?0:maximum([depth(i) for i in c.container])
depth(node::DomainDecompositionNode) = max(depth(node.Mid), depth(node.Other))+1

"Interpret i in a binary way if i=011. Take the route Other, Other, Mid to reach the required node"
get_nodes(node::DomainDecompositionNode, i::Int=depth(node)) =  get_nodes(iseven(i) ? node.Mid : node.Other, i>>1)
get_nodes(c::DomainDecompositionNodeContainer, i::Int=depth(c)) = vcat([get_nodes(ci, i) for ci in c.container]...)
get_nodes(leaf::DomainDecompositionLeaf, i::Int=depth(leaf)) = leaf

"All leave container below the given node"
leave_containers(c::DomainDecompositionNodeContainer) = is_leaf_container(c) ? c :  vcat([leave_containers(i) for i in c.container]...)
leave_containers(node::DomainDecompositionNode) = vcat(leave_containers(node.Mid), leave_containers(node.Other) )

is_leaf_container(c::DomainDecompositionNodeContainer) = eltype(c.container) <: DomainDecompositionLeaf

"""
All nodes that are below `no_other_turns` other turns

Example:
```jldoctest
julia> get_all_nodes(tree, 2, 3)
```
returns the nodes which are obtained by
```jldoctest
julia> get_nodes(tree, 3, 3)
julia> get_nodes(tree, 5, 3)
julia> get_nodes(tree, 6, 3)
```
where 3, 5, 6 are number os 3 bit of which 2 are 1.
"""
function get_all_nodes(tree, no_other_turns::Int, depth::Int=depth(tree))
    r = DomainDecompositionLeaf[]
    for i in get_combinations(no_other_turns, depth)
        push!(r, get_nodes(tree, i)...)
    end
    r
end

# required for pretty printing.
"Get the children of a node"
children(::DomainDecompositionLeaf) = ()
children(node::DomainDecompositionNode) = (node.Mid, node.Other)
children(c::DomainDecompositionNodeContainer) = c.container

#  required for pretty printing.
printnode(io::IO, d::DomainDecompositionNodeContainer) = Base.print_with_color(:blue, io, "C")
printnode(io::IO, d::DomainDecompositionNode) = print(io, "node")
printnode(io::IO, f::FrameFun.DomainDecompositionLeaf) = print(io, "grid: "*string(size(FrameFun.DMZ(f))))

import Base:split
# Split the leaf by using its parameters
split(leaf::DomainDecompositionLeaf, basis, domain; options...) = split(DMZ(leaf), BasisFunctions.grid(leaf), basis, domain; options...)
split(DMZ::BasisFunctions.AbstractGrid, basis, domain; options...) = split(DMZ, nothing, basis, domain; options...)
function Base.split(DMZ::MaskedGrid, grid::Union{Void,MaskedGrid}, basis, domain; splitting=true, options...)
    # Z is a grid which intersects domain and which has the width of the span of a basis element
    if length(DMZ) < 1000
        return FrameFun.DomainDecompositionNode(FrameFun.DomainDecompositionNodeContainer(DMZ, grid), FrameFun.DomainDecompositionNodeContainer(FrameFun.DomainDecompositionLeaf[]))
    end
    split_DMZ = FrameFun.subgrid(DMZ, domain)
    Z = FrameFun.boundary_support_grid(basis, split_DMZ, DMZ)

    # If Z has no elements, no decompsition is possible
    if length(Z) == 0
        #  A node in which the Other part is empty and the Mid part is the input.
        return FrameFun.DomainDecompositionNode(FrameFun.DomainDecompositionNodeContainer(DMZ, grid), FrameFun.DomainDecompositionNodeContainer(FrameFun.DomainDecompositionLeaf[]))
    end
    # Split the left by Z
    # if splitting
        grids = split(DMZ - Z)
    # else
        # grids = [DMZ-Z]
    # end

    # Determine the size of the region of influence of spines overlapping with grids
    DMZs = map(x->FrameFun.boundary_support_grid(basis, x, DMZ), grids)
    # Embed Z into split_DMZ_DMZ to ensure that coefficients (not only evaluations) are solved accurately in Z
    split_DMZ_DMZ = FrameFun.boundary_support_grid(basis, Z, DMZ)
    if length(split_DMZ_DMZ) == length(DMZ)
        #  A node in which the Other part is empty and the Mid part is the input.
        return FrameFun.DomainDecompositionNode(FrameFun.DomainDecompositionNodeContainer(DMZ, grid), FrameFun.DomainDecompositionNodeContainer(FrameFun.DomainDecompositionLeaf[]))
    end


    # It might be that also split_DMZ_DMZ can be split into different domains
    # if splitting
        split_DMZ_DMZs = split(split_DMZ_DMZ)
    # else
        # split_DMZ_DMZs = [split_DMZ_DMZ]
    # end

    # If grid is nothing we are splitting a domain in the path Mid, ..., Mid and there all splines can be used to solve.
    if grid==nothing
        FrameFun.DomainDecompositionNode(FrameFun.DomainDecompositionNodeContainer(split_DMZ_DMZs),
            FrameFun.DomainDecompositionNodeContainer(DMZs, grids))
    else
        # grid indicates which splines can be used to solve the input problem.
        # Make sure that the splines that can be used in the splitted problem are the same
        m = BasisFunctions.mask(grid)
        split_gs = map(x->MaskedGrid(supergrid(DMZ), BasisFunctions.mask(x) .& m), split_DMZ_DMZs)
        gs = map(x->MaskedGrid(supergrid(DMZ), BasisFunctions.mask(x) .& m), grids)
        FrameFun.DomainDecompositionNode(FrameFun.DomainDecompositionNodeContainer(split_DMZ_DMZs, split_gs),
            FrameFun.DomainDecompositionNodeContainer(DMZs, gs))
    end
end

function split(node::DomainDecompositionNode, basis, domain; options...)
    other = [split(g, basis, domain; options...) for g in node.Other.container]
    mid = [split(g, basis, domain; options...) for g in node.Mid.container]
    DomainDecompositionNode(DomainDecompositionNodeContainer(mid), DomainDecompositionNodeContainer(other))
end

"""
Decompose the grids in the leaves of the tree using the domain.
"""
function split!(tree, basis, domain; options...)
    for c in leave_containers(tree)
        r = Array{DomainDecompositionNode}(length(c.container))
        for i in 1:length(c.container)
            r[i] = split(c.container[i], basis, domain; options...)
        end
        c.container =  r
    end
    tree
end

# Util functions for plotting domain decomposition trees
using RecipesBase
@recipe function f(leaf::DomainDecompositionLeaf, level=4)
    legend --> false
    markersize --> 6-level
    DMZ(leaf)
end
@recipe function f(c::DomainDecompositionNodeContainer, level=4)
    legend --> false
    for i in c.container
        @series begin
            i, level
        end
    end
end
@recipe function f(node::DomainDecompositionNode, level=4)
    legend --> false
    @series begin
        node.Mid, level-2
    end
    @series begin
        node.Other, level
    end
end

# return the coordinates of places interesting for decomposing a domain
# Used in split_domain_Nd
function domain_grid_Nd(basis::Dictionary, gamma::AbstractGrid, ranges::Domain; options...)
    bb = boundingbox(ranges)
    r = [(ai, bi) for (ai, bi) in zip(leftendpoint(bb), rightendpoint(bb))]
    [domain_grid_1d(d, basis, gamma, r; options...) for d in 1:dimension(basis)]
end

"The coordinates of places interesting for decomposing a domain in 1 dimension"
function domain_grid_1d(d::Int, basis::Dictionary, gamma::AbstractGrid, ranges;
        shift::Bool=false, factor::Real=3, side_length=1000^(1/(dimension(basis)-1)), options...)
    basis_coarseness = BasisFunctions.support_length_of_compact_function(element(basis, d))
    grid_coarsness = BasisFunctions.stepsize(element(gamma, d))
    step = max(factor*basis_coarseness, grid_coarsness*side_length)
    left = min(mean(ranges[d]),minimum(ranges[d]) + step)
    right = max(mean(ranges[d]),maximum(ranges[d]) - step)
    N = floor(Int, max(0,(right-left)/(step)))
    if !shift
        if N == 0 || N == 1
            return tuple((left+right)/2)
        end
        return tuple(linspace(left,right,N)...)
    else
        if N == 0 || N==1
            error("shift not possible")
        end
        if N == 2
            return tuple((left+right)/2)
        end
        return tuple(linspace(left+step, right-step,N-1)...)
    end
end


"""
A domain in `gamma`, that splits `gamma` in dimenions `dim` at the coordinate 'mid'
"""
function split_domain_1d(gamma::AbstractGrid, dim::Int, mid::ELT) where {ELT <: Real}
    a = leftendpoint(gamma); b = rightendpoint(gamma)
    dx = BasisFunctions.stepsize(elements(gamma)[dim])
    D = length(dx)
    if dim ==1
        split_domain = interval(mid-dx/2,mid+dx/2)×Domains.ProductDomain([interval(a[i],b[i]) for i in 2:length(a)]...)
    elseif dim==D
        split_domain = ProductDomain([interval(a[i],b[i]) for i in 1:length(a)-1]...)×interval(mid-dx/2,mid+dx/2)
    else
        split_domain = ProductDomain([interval(a[i],b[i]) for i in 1:dim-1]...)×interval(mid-dx/2,mid+dx/2)×ProductDomain([interval(a[i],b[i]) for i in dim+1:length(a)]...)
    end
    split_domain
end

"""
A domain in `gamma` that splits `gamma` in dimenions `dim` at the coordinates 'mid'.
"""
split_domain_1d(gamma::AbstractGrid, dim::Int, mid::NTuple{N,ELT}) where {N,ELT<:Real} =
    UnionDomain([split_domain_1d(gamma, dim, mid[d]) for d in 1:length(mid)]...)

"""
A domain in `gamma` that splits `gamma` in  coordinates in 'mid[d]' for dimension `d`.
"""
split_domain_Nd(gamma::AbstractGrid, mids::Vector) =
    UnionDomain([split_domain_1d(gamma, d, m) for (d,m) in enumerate(mids)]...)

"""
A domain in `gamma` that splits `gamma` in the equispaced inside the ranges
"""
split_domain_Nd(basis::Dictionary, gamma::AbstractGrid, ranges::Domain; options...) =
    split_domain_Nd(gamma, domain_grid_Nd(basis, gamma, ranges; options...))

"""
A domain decomposition tree for the grid 'DMZ' using the coarseness of elements in 'basis'
and the coarseness of `gamma`. Splitting occurs in the bounding box of `ranges`
"""
function create_tree(DMZ::AbstractGrid, basis::Dictionary, gamma::AbstractGrid, ranges::Domain;
        depth::Int=dimension(gamma), options...)
    split_domain =FrameFun.split_domain_Nd(basis, gamma, ranges; options...)
    # println(split_domain)
    D = dimension(gamma)
    domain_array = elements(split_domain)
    @assert length(domain_array) == D
    tree = split(DMZ, basis, domain_array[1]; options...)
    for d in 2:min(D, depth)
        if d == 2
            tree = split(tree, basis, domain_array[d]; options...)
        else
            FrameFun.split!(tree, basis, domain_array[d]; options...)
        end
    end
    tree
end


create_tree(basis::Dictionary, gamma::AbstractGrid, omega::AbstractGrid, domain::Domain; options...) =
    create_tree(boundary_support_grid(basis, boundary_grid(gamma, domain), omega), basis, gamma, domain; options...)

"All integers of 'depth' bits that have 'no' ones"
function get_combinations(no::Int, depth::Int)
    a = zeros(Int, depth)
    a[1:no] = 1
    r = Set{Int}()
    for p in permutations(a)
        push!(r, _int(p))
    end
    collect(r)
end

"Convert array of ones and zeros to integer. [1,1,0] => 110b = 6. "
function _int(a::Vector{Int})
    i = 0
    for ai in a
        i = 2*i+ai
    end
    i
end

"""
Solving the problem Ax=b where b=f(x_i) and x_i∈G.
The tree contains a domain decomposition for the boundary of the grid G.
"""
function solve(b::Vector, A::AbstractOperator, tree::AbstractDomainDecompositionNode;
        options...)
    omega = BasisFunctions.grid(dest(A))
    gamma = supergrid(omega)
    basis = superdict(dictionary(src(A)))
    solve(b, tree, basis, gamma, omega; options...)
end

function solve(b::Vector, A::AbstractOperator, tree::AbstractDomainDecompositionNode,
        basis::Dictionary, gamma::AbstractGrid, omega::AbstractGrid; options...)
    x0 = zeros(basis)
    b_ext = restriction_operator(gridspace(gamma), gridspace(omega))'*b
    s = GridSamplingOperator(gridspace(gamma))
    solve!(x0, A, b_ext, tree, basis, gamma, omega, s*basis, depth(tree), copy(b_ext), zeros(length(b)))
end

function solve!(x0::Array{ELT}, A::AbstractOperator, b_ext::Array{ELT},
        tree::FrameFun.AbstractDomainDecompositionNode, basis::Dictionary,
        gamma::AbstractGrid, omega::AbstractGrid, OP::AbstractOperator,
        depth::Int, b_ext_copy::Array{ELT}, bnew::Vector{ELT};
        cutoff=FrameFun.default_cutoff(OP)) where {ELT<:Number}
    to = 2^depth-1
    for i in 0:to
        n = 0
        for leaf in FrameFun.get_nodes(tree, i)
            GR = restriction_operator(gridspace(gamma), gridspace(FrameFun.DMZ(leaf)))
            if FrameFun.grid(leaf) == nothing
                boundary_indices = FrameFun.boundary_element_indices(basis, FrameFun.DMZ(leaf))
            else
                boundary_indices = FrameFun.boundary_element_indices(basis, FrameFun.grid(leaf))
            end
            FR = IndexRestrictionOperator(Span(basis), Span(basis[boundary_indices]), boundary_indices)
            a = matrix(GR*OP*FR')
            # println(size(a))
            y0 = LAPACK.gelsy!(a, GR*b_ext, cutoff)[1]
            # println( norm(GR*OP*FR'*y0- GR*b_ext ))
            x0[FR.subindices] .+= y0
            n += 1
        end
        if i < to && n > 0
            apply!(A, bnew, x0)
            for (i_i,i) in enumerate(subindices(omega))
                b_ext[i] = b_ext_copy[i] - bnew[i_i]
            end
        end
    end
    x0
end
