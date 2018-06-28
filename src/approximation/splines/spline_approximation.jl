SYSTEM_SIZE=2000

"""
Number of basis elements overlapping with a point.
"""
no_overlapping_elements(dict::Dictionary) = ceil(Int,BasisFunctions.support_length_of_compact_function(dict)*length(dict))

"""
Assign a sequence number to each element in bin.

Elements with all numbers even are to be solved first and are assigned 1.
Elements with all numbers odd are to be solved last.
"""
function assign_sequence_nro(bins)
    L = length(bins)
    M = length(bins[1])
    seq = Array{Int}(L)
    a = Array{Bool}(M)
    b = Array{Int}(M)
    for i in 1:L
        b .= bins[i]
        a .= iseven.(b)
        seq[i] = M+1-sum(a)
    end
    seq
end

"""
Clasifies the coefficient indices activated in `coefficient_mask` in `depth` dimensions.
"""
function classified_indices(coefficient_mask::AbstractArray{Bool}, primal::TensorProductDict, gamma::ProductGrid, depth::Int)
    sprimal = size(primal)
    # The cartesian indices of the activated coefficients
    cart_indices = [CartesianIndex(ind2sub(sprimal, i) ) for i in find(coefficient_mask)]
    cart_indices_matrix = zeros(Int,length(cart_indices), depth)
    for i in 1:depth
        # Transfrom one dimensian of the cartesian indices to an array
        cart_indices_matrix[:,i] .= getindex.(cart_indices,i)
        # Determin the space between two coefficient regions in the ith dimension
        # Take into account the spacing of the collocation points and the width
        # or the primal and the dual basis.
        # (it scales, but other methods may be better)
        g1d = element(gamma,i)
        primal1d = element(primal, i)
        dual1d = BasisFunctions.wavelet_dual(primal1d)
        OS = cld(length(g1d),length(primal1d))
        no_samples_in_1d = cld(SYSTEM_SIZE,OS*no_overlapping_elements(dual1d))
        no_coeffs_in_1d_other = ceil(Int, fld(no_samples_in_1d, OS)^(1/(max(1,depth-1))))
        no_coeffs_in_1d_mid = ceil(Int, no_overlapping_elements(primal1d)^(1/(max(1,depth-1))))

        # Divide all coefficients in a single dimension in an even n.o. partitions (m).
        m = cld(length(primal1d),max(no_coeffs_in_1d_mid, no_coeffs_in_1d_other))
        m = isodd(m) ? m+1 : m
        # +1 to be on the save side.
        m = (length(primal1d)+1)/m
        # Reuse the array to minimize allocation
        cart_indices_matrix[:,i] .= Int.(cld.(cart_indices_matrix[:,i], m))
    end
    # Create a vector instead of a matrix.
    # Then reducing dimensions is not necessary in methods using this output.
    cart_indices, [tuple(cart_indices_matrix[i,:]...) for i in 1:size(cart_indices,1)]
end

function azselection_restriction_operators(fplatform::BasisFunctions.GenericPlatform, i; options...)
    platform = fplatform.super_platform
    s = sampler(fplatform, i)
    omega = grid(s)
    gamma = supergrid(omega)
    domain = FrameFun.domain(primal(fplatform, i))
    azselection_restriction_operators(primal(platform, i), gamma, omega, domain)
end

azselection_restriction_operators(primal::Union{BasisFunctions.WaveletBasis,BasisFunctions.WaveletTensorDict}, gamma::AbstractGrid, omega::AbstractGrid, domain::Domains.Domain) =
    azselection_restriction_operators(primal, BasisFunctions.wavelet_dual(primal), gamma, omega, domain)

azselection_restriction_operators(primal::Dictionary, gamma::AbstractGrid, omega::AbstractGrid, domain::Domains.Domain) =
    azselection_restriction_operators(primal, primal, gamma, omega, domain)


"""
Grid restriction and dictionary restriction operators to restrict the system when approximating with compact frame elements.

`primal` and `dual` should be the not restricted basis, `gamma` the (oversampled) grid of `primal`, `omega` is `gamma` restricted to `domain`
"""
function azselection_restriction_operators(primal::Dictionary, dual::Dictionary, gamma::AbstractGrid, omega::MaskedGrid, domain::Domains.Domain)
    bound = FrameFun.boundary_grid(gamma, domain)
    coefficient_mask = BasisFunctions.coefficient_index_mask_of_overlapping_elements(dual, bound)
    _azselection_restriction_operators(primal, gamma, omega, coefficient_mask)
end

function _azselection_restriction_operators(primal::Dictionary, gamma::AbstractGrid, omega::AbstractGrid, coefficient_mask)
    grid_mask = BasisFunctions.grid_index_mask_in_element_support(primal, gamma, coefficient_mask)
    grid_mask .= grid_mask .& FrameFun.mask(omega)
    __azselection_restriction_operators(primal, gamma, omega, grid_mask)
end

function __azselection_restriction_operators(primal::Dictionary, gamma::AbstractGrid, omega::AbstractGrid, grid_mask)
    DMZ = gamma[grid_mask]
    system_coefficient_mask = BasisFunctions.coefficient_index_mask_of_overlapping_elements(primal, DMZ)
    gr = restriction_operator(gridbasis(omega, coeftype(primal)), gridbasis(DMZ, coeftype(primal)))
    dr = restriction_operator(primal, system_coefficient_mask)
    dr, gr
end

"""
A grid that contains the points of `omega_grid` that are not evaluated to zero by the elements that overlap with boundary_grid.
"""
boundary_support_grid(basis, boundary_grid::Union{MaskedGrid,IndexSubGrid}, omega_grid::Union{MaskedGrid,IndexSubGrid}) =
    boundary_support_grid(basis, basis, boundary_grid, omega_grid)

function boundary_support_grid(basis::Dictionary, dual::Dictionary, boundary_grid::Union{MaskedGrid,IndexSubGrid}, omega_grid::Union{MaskedGrid,IndexSubGrid})
    boundary_element_m = BasisFunctions.coefficient_index_mask_of_overlapping_elements(dual, boundary_grid)
    gamma = supergrid(omega_grid)
    m = BasisFunctions.grid_index_mask_in_element_support(basis, gamma, boundary_element_m)
    m .= m .& FrameFun.mask(omega_grid)
    MaskedGrid(gamma,m)
end

spline_util_restriction_operators(platform::BasisFunctions.GenericPlatform, i) =
    spline_util_restriction_operators(primal(platform, i), sampler(platform, i))

spline_util_restriction_operators(dict::Dictionary, sampler::GridSamplingOperator) =
    spline_util_restriction_operators(dict, grid(sampler))

spline_util_restriction_operators(frame::ExtensionFrame, grid::Union{IndexSubGrid,MaskedGrid}) =
    spline_util_restriction_operators(superdict(frame), grid, supergrid(grid), Domains.domain(frame))

spline_util_restriction_operators(dict::TensorProductDict{N,NTuple{N1,DT},S,T}, grid::AbstractGrid) where {N,N1,DT<:ExtensionFrame,S,T} =
    spline_util_restriction_operators(FrameFun.flatten(dict), grid)

spline_util_restriction_operators(dict::Dictionary, omega::AbstractGrid, grid::AbstractGrid, domain::Domain) =
    spline_util_restriction_operators(dict, omega, boundary_grid(grid, domain))

"""
Frame restriction operator and grid restriction operator.
The former restricts `dict` to the elements that overlap with the boundary and
the latter restricts `omega` to the points in the span of the dict elements that
overlap with the boundary.
"""
spline_util_restriction_operators(dict::Dictionary, omega::AbstractGrid, boundary::AbstractGrid) =
    _spline_util_restriction_operators(dict, omega, boundary_support_grid(dict, boundary, omega))

function BasisFunctions.grid_restriction_operator(src::Dictionary, dest::Dictionary, src_grid::Union{IndexSubGrid,MaskedGrid}, dest_grid::MaskedGrid)
    @assert supergrid(src_grid) == supergrid(dest_grid)
    IndexRestrictionOperator(src, dest, FrameFun.relative_indices(dest_grid,src_grid))
end

function _spline_util_restriction_operators(dict::Dictionary, grid::AbstractGrid, DMZ::AbstractGrid)
    boundary_indices = BasisFunctions.coefficient_indices_of_overlapping_elements(dict, DMZ)
    frame_restriction = IndexRestrictionOperator(dict, dict[boundary_indices], boundary_indices)
    grid_restriction = BasisFunctions.restriction_operator(gridbasis(grid), gridbasis(DMZ))
    frame_restriction, grid_restriction
end

"""
Frame restriction operator and grid restriction operator.

The former restricts `dict` to the elements that overlap with the boundary and
the latter restricts `grid` to the points in the span of the dict elements that
overlap with a region defined as the span of the dict elements that overlap with the boundary.
"""
function _spline_util_restriction_operators(dict::Dictionary, grid::AbstractGrid, boundary::AbstractGrid, DMZ::AbstractGrid)
    boundary_indices = BasisFunctions.coefficient_indices_of_overlapping_elements(dict, boundary)
    frame_restriction = IndexRestrictionOperator(dict, dict[boundary_indices], boundary_indices)
    grid_restriction = BasisFunctions.restriction_operator(gridbasis(grid), gridbasis(DMZ))
    frame_restriction, grid_restriction
end

estimate_plunge_rank(src::BSplineTranslatesBasis, domain::Domain, dest::GridBasis) =
    estimate_plunge_rank(src, domain, grid(dest))
estimate_plunge_rank(src::BSplineTranslatesBasis, domain::Domain, grid::Union{IndexSubGrid,MaskedGrid}) =
    estimate_plunge_rank_spline(src, domain, supergrid(grid))

estimate_plunge_rank(src::TensorProductDict{N,DT,S,T}, domain::Domain, dest::GridBasis) where {N,DT<:NTuple{N1,BasisFunctions.BSplineTranslatesBasis} where {N1},S,T} =
    estimate_plunge_rank(src, domain, grid(dest))

estimate_plunge_rank(src::TensorProductDict{N,DT,S,T}, domain::Domain, grid::MaskedGrid) where {N,DT<:NTuple{N1,BasisFunctions.BSplineTranslatesBasis} where {N1},S,T} =
    estimate_plunge_rank_spline(src, domain, supergrid(grid))

function estimate_plunge_rank_spline(src, domain::Domain, grid::AbstractGrid)
    boundary = boundary_grid(grid, domain)
    sum(BasisFunctions.coefficient_index_mask_of_overlapping_elements(src, boundary))
end

# The grid on the boundary of omega
divide_and_conquer_restriction_operators(omega::AbstractGrid, gamma::AbstractGrid,
        basis::Dictionary, domain::Domains.Domain, dim::Int, range::AbstractVector) =
    divide_and_conquer_restriction_operators(omega::AbstractGrid, gamma::AbstractGrid,
        basis::Dictionary, FrameFun.boundary_support_grid(basis, boundary_grid(gamma, domain), omega)::AbstractGrid, dim::Int, range::AbstractVector)

function intersect_boundary(DMZ::AbstractGrid, gamma::AbstractGrid, dim::Int, range)
    a = leftendpoint(gamma); b = rightendpoint(gamma)
    dx = BasisFunctions.stepsize(elements(gamma)[dim])
    mid = mean(range)
    D = dimension(gamma)
    # Some domain restriction one dimension of the bounding box
    if dim ==1
        split_domain = interval(mid-dx/2,mid+dx/2)×Domains.ProductDomain([interval(a[i],b[i]) for i in 2:D]...)
    elseif dim==D
        split_domain = ProductDomain([interval(a[i],b[i]) for i in 1:D-1]...)×interval(mid-dx/2,mid+dx/2)
    else
        split_domain = ProductDomain([interval(a[i],b[i]) for i in 1:dim-1]...)×interval(mid-dx/2,mid+dx/2)×ProductDomain([interval(a[i],b[i]) for i in dim+1:D]...)
    end

    # The grid on the intersection of split_domain and the boundary
    split_grid = FrameFun.subgrid(DMZ, split_domain)
    if length(split_grid) == 0
        error("check your ranges for splitting the domain.")
    end
    split_grid
end

function FrameFun.divide_and_conquer_restriction_operators(omega::AbstractGrid,
        gamma::AbstractGrid, basis::Dictionary, DMZ::AbstractGrid, dim::Int, range::AbstractVector)
    split_grid = FrameFun.intersect_boundary(DMZ, gamma, dim, range)
    #  The points on the support of the functions overlapping with split_grid
    mid = FrameFun.boundary_support_grid(basis, split_grid, DMZ)
    # Look at the grid away from the intersection and split in two disconnected parts.
    left_over = DMZ-mid
    if length(left_over) != 0
        left, right = split(left_over)
        leftDMZ = FrameFun.boundary_support_grid(basis, left, DMZ)
        rightDMZ = FrameFun.boundary_support_grid(basis, right, DMZ)

        if !(.5 < length(left)/length(right) < 2)
            warn("ratio between the separated grids is not distributed evenly: $(length(left)) and $(length(right))")
        end

        BR0, DMZ_R0 = FrameFun._spline_util_restriction_operators(basis, omega, FrameFun.boundary_support_grid(basis, mid, DMZ))
        BR1, DMZ_R1 = FrameFun._spline_util_restriction_operators(basis, omega, left, leftDMZ)
        BR2, DMZ_R2 = FrameFun._spline_util_restriction_operators(basis, omega, right, rightDMZ)

        BR0, BR1, BR2, DMZ_R0, DMZ_R1, DMZ_R2
    else
        warn("No grid splitting possible")
        BR, DMZ_R = FrameFun._spline_util_restriction_operators(basis, omega, DMZ)
    end
end


function split_DMZ(DMZ::AbstractGrid, basis::Dictionary, gamma::AbstractGrid, ranges::Domain; options...)
    split_grid = intersect_DMZ(DMZ::AbstractGrid, basis, gamma::AbstractGrid, ranges; options...)
    split_grid_DMZ = FrameFun.boundary_support_grid(basis, split_grid, DMZ)
    left_over = DMZ - split_grid_DMZ
    if length(left_over) == 0
        return [DMZ], AbstractGrid[]
    end
    FrameFun.split_in_IndexGrids(FrameFun.boundary_support_grid(basis, split_grid_DMZ, DMZ)), FrameFun.split_in_IndexGrids(left_over)
end


function intersect_DMZ(DMZ::AbstractGrid, basis::Dictionary, gamma::AbstractGrid, ranges::Domain; verbose=false, options...)
    domain_grid = domain_grid_Nd(basis, gamma, ranges; options...)
    split_domain = split_domain_Nd(gamma, domain_grid)
    verbose && println(split_domain)

    # The grid on the intersection of split_domain and the boundary
    split_grid = FrameFun.subgrid(DMZ, split_domain)
    if length(split_grid) == 0
        error("check your mids for splitting the domain.")
    end
    split_grid
end

function divide_and_conquer_N_util_operators(fplatform::BasisFunctions.Platform, i; ranges=nothing, options...)
    platform = fplatform.super_platform
    basis = primal(platform, i)
    S = sampler(fplatform, i)
    s = sampler(platform, i)
    domain = FrameFun.domain(primal(fplatform, i))
    gamma = BasisFunctions.grid(s)
    omega = BasisFunctions.grid(S)
    (ranges==nothing) && (ranges=domain)
    FrameFun.divide_and_conquer_N_util_operators(omega, gamma, basis, domain, ranges; options...)
end


# The grid on the boundary of omega
divide_and_conquer_N_util_operators(omega::AbstractGrid, gamma::AbstractGrid, basis::Dictionary, domain::Domains.Domain, ranges::Domain; options...) =
    divide_and_conquer_N_util_operators(omega::AbstractGrid, gamma::AbstractGrid, basis::Dictionary, FrameFun.boundary_support_grid(basis, boundary_grid(gamma, domain), omega)::AbstractGrid, ranges::Domain; options...)

function FrameFun.divide_and_conquer_N_util_operators(omega::AbstractGrid, gamma::AbstractGrid, basis::Dictionary, DMZ::AbstractGrid, ranges::Domain; recur=nothing, options...)
    DMZs, GRs = FrameFun.split_DMZ(DMZ, basis, gamma, ranges; options...)
    OP = GridSamplingOperator(gridbasis(gamma))*basis
    ops = []
    (recur == nothing) && (recur = dimension(basis)==3 ? 2 :1 )
    if recur >= 1
        ops = push!(ops, FrameFun.util_operators(OP, DMZ, DMZs, GRs, gamma, basis)...)
    end
    if recur >= 2
        @assert length(DMZs)==1
        try
        SplitDMZs, SplitGRs = FrameFun.split_DMZ(DMZs[1], basis,  gamma, ranges; shift=true)
        ops = push!(ops, FrameFun.util_operators(OP, DMZs[1], SplitDMZs, SplitGRs, gamma, basis)...)
        catch
            warn("no extra splitting was possible")
        end
    end
    tuple(ops...)
end

function util_operators(OP, DMZ, DMZs, GRs, gamma, basis)
    A0 = Array{CompositeOperator}(length(DMZs))
    GR0 = Array{IndexRestrictionOperator}(length(DMZs))
    for (i,dmz) in enumerate(DMZs)
        frame_restriction, grid_restriction = FrameFun._spline_util_restriction_operators(basis, gamma, dmz)
        A0[i] = grid_restriction*OP*frame_restriction'
        GR0[i] = grid_restriction
    end

    A1 = Vector{CompositeOperator}(length(GRs))
    GR1 = Vector{IndexRestrictionOperator}(length(GRs))
    for (i,gr) in enumerate(GRs)
        dmz = FrameFun.boundary_support_grid(basis, gr, DMZ)
        frame_restriction, grid_restriction = FrameFun._spline_util_restriction_operators(basis, gamma, gr, dmz)
        A1[i] = grid_restriction*OP*frame_restriction'
        GR1[i] = grid_restriction
    end
    A0, GR0, A1, GR1
end
