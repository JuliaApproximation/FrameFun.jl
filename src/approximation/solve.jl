

# Function with equal functionality, but allocating memory
function az_solve(b, A::DictionaryOperator, Zt::DictionaryOperator, util...; use_plunge=true,
        cutoff = default_cutoff(A), trunc = truncatedsvd_solve, verbose=false, afirst=true, options...)
    if afirst
        if use_plunge
            P = plunge_operator(A, Zt)
            x2 = trunc(P*b, P*A, util...; cutoff=cutoff, verbose=verbose, options...)
        else
            x2 = trunc(b, A, util... ; cutoff=cutoff, verbose=verbose, options...)
        end
        x1 = Zt*(b-A*x2)
    else
        x1 = Zt*b
        x2 = trunc(b-A*x1, A, util... ; cutoff=cutoff, verbose=verbose, options...)
    end
    x1 .= x1 .+ x2
    x1
end

function az_solve(platform::BasisFunctions.Platform, i, f::Function; R=0, options...)
    a = A(platform, i)
    zt = Zt(platform, i)
    s = sampler(platform, i)
    (R == 0) && (R=estimate_plunge_rank(a))
    az_solve(s*f, a, zt; R=R, options...)
end

# Function with equal functionality, but allocating memory
azs_solve(b, A::DictionaryOperator, Zt::DictionaryOperator, RD::DictionaryOperator, SB::DictionaryOperator;
        trunc = restriction_solve, use_plunge=false, info=false, options...) =
    az_solve(b, A, Zt, RD, SB; trunc=trunc, use_plunge=use_plunge, options...)

function azs_solve(fplatform::BasisFunctions.Platform, i, f::Function; info=false,options...)
    a = A(fplatform, i; options...)
    zt = Zt(fplatform, i; options...)
    platform = fplatform.super_platform
    s = sampler(fplatform, i)
    omega = grid(s)
    gamma = supergrid(omega)
    domain = FrameFun.domain(src(a))
    frame_restriction, grid_restriction = azselection_restriction_operators(primal(platform, i), gamma, omega, domain)
    if info
        restriction_info(s*f, a, frame_restriction', grid_restriction)
    else
        azs_solve(s*f, a, zt, frame_restriction', grid_restriction; options...)
    end
end

function az_decomposition_solve(fplatform::BasisFunctions.Platform, i, f::Function;
        depth=nothing, info=false, fig=false, no_blocks=nothing, options...)
    platform = fplatform.super_platform
    a = A(fplatform, i)
    zt = Zt(fplatform, i)
    S = sampler(fplatform, i)

    # The grid on Gamma
    gamma = grid(sampler(platform, i))
    # The grid on Omega
    omega = grid(S)

    dom = domain(primal(fplatform, i))
    basis = primal(platform, i)
    (depth==nothing) && (depth=dimension(basis))
    dual = BasisFunctions.wavelet_dual(basis)
    bound = FrameFun.boundary_grid(gamma, dom)
    boundary_coefficient_mask = BasisFunctions.coefficient_index_mask_of_overlapping_elements(dual, bound)
    cart_indices, c_indices = classified_indices(boundary_coefficient_mask, basis, gamma, depth; no_blocks=no_blocks)

    if info
        decomposition_info(S*f, a, cart_indices, c_indices)
    elseif fig
        decomposition_plot(a, cart_indices, c_indices)
    else
        az_solve(S*f, a, zt, cart_indices, c_indices; trunc=decomposition_solve, use_plunge=false, options...)
    end
end

function timed_az_decomposition_solve(fplatform::BasisFunctions.Platform, i, f::Function;
        afirst=false, depth=nothing, no_blocks=nothing, verbose=false, options...)
    platform = fplatform.super_platform
    t1 = @timed begin
        platform = fplatform.super_platform
        a = A(fplatform, i)
        zt = Zt(fplatform, i)
        S = sampler(fplatform, i)

        # The grid on Gamma
        gamma = grid(sampler(platform, i))
        # The grid on Omega
        omega = grid(S)

        dom = domain(primal(fplatform, i))
        basis = primal(platform, i)
        (depth==nothing) && (depth=dimension(basis))
        dual = BasisFunctions.wavelet_dual(basis)
        b = S*f
    end
    t2 = @timed bound = FrameFun.boundary_grid(gamma, dom)
    t3 = @timed boundary_coefficient_mask = BasisFunctions.coefficient_index_mask_of_overlapping_elements(dual, bound)
    t4 = @timed cart_indices, c_indices = classified_indices(boundary_coefficient_mask, basis, gamma, depth; no_blocks=no_blocks)
    if afirst
        # t5 = @timed x2 = decomposition_solve(b, a, cart_indices, c_indices ; verbose=verbose, options...)
        # t6 = @timed nothing
        t5 = @timed matrices = decomposition_matrices(a, cart_indices, c_indices ; verbose=verbose, options...)
        t6 = @timed x2 = decomposition_solve_matrices(b,  a, cart_indices, c_indices, matrices ; verbose=verbose, options...)
        t7 = @timed x1 = zt*(b-a*x2)
    else
        t7 = @timed x1 = zt*b
        t5 = @timed matrices = decomposition_matrices(a, cart_indices, c_indices ; verbose=verbose, options...)
        t6 = @timed x2 = decomposition_solve_matrices(b-a*x1,  a, cart_indices, c_indices, matrices ; verbose=verbose, options...)
    end
    t8 = @timed x1 .= x1 .+ x2
    info("error is $(norm(a*x1-b))")
    x1, [t1,t2,t3,t4,t5,t6,t7,t8]
end

function decomposition_matrices(A::DictionaryOperator{ELT}, cart_indices, classified_indices; cutoff=default_cutoff(A), verbose=false, info=false, options...) where {ELT}
    bins = unique(classified_indices)
    # assign sequence numbers to each of the bins.
    seq_nr = assign_sequence_nro(bins)

    primal = basis(src(A))
    omega = grid(dest(A))
    g = supergrid(omega)
    x1 = zeros(primal)
    j = 1
    r = Array{Matrix{ELT}}(length(bins))
    for d in 1:1<<length(bins[1])
        # first solve all parts with a low sequence number
        bpart = bins[find(seq_nr.==d)]
        # Each of the parts with an equal sequence number can be solved independently
        verbose && println("$(length(bpart)) parts in $(d)th flow")
        for i in 1:size(bpart,1)
            mo = cart_indices[find(classified_indices.==bpart[i,:])]
            xx, yy = FrameFun._azselection_restriction_operators(primal, g, omega, mo)
            verbose && println("\t$(i)\t has size ($(size(yy,1)),$(size(xx,1)))")
            op = yy*A*xx'
            a = matrix(op)
            r[j] = a
            j += 1
        end
    end
    r
end


function decomposition_solve_matrices(b::Vector, A::DictionaryOperator{ELT}, cart_indices, classified_indices, matrices; cutoff=default_cutoff(A), verbose=false, info=false, options...) where {ELT}
    bins = unique(classified_indices)
    # assign sequence numbers to each of the bins.
    seq_nr = assign_sequence_nro(bins)

    primal = basis(src(A))
    omega = grid(dest(A))
    g = supergrid(omega)
    x1 = zeros(primal)
    b_ =  copy(b)
    t = similar(x1)

    j = 1
    for d in 1:1<<length(bins[1])
        # first solve all parts with a low sequence number
        bpart = bins[find(seq_nr.==d)]
        # Each of the parts with an equal sequence number can be solved independently
        verbose && println("$(length(bpart)) parts in $(d)th flow")
        for i in 1:size(bpart,1)
            mo = cart_indices[find(classified_indices.==bpart[i,:])]
            xx, yy = FrameFun._azselection_restriction_operators(primal, g, omega, mo)
            verbose && println("\t$(i)\t has size ($(size(yy,1)),$(size(xx,1)))")
            a = matrices[j]
            j += 1
            y = LAPACK.gelsy!(a, yy*b_, cutoff)[1]
            # x1 = x1 + xx'*y
            apply!(xx', t, y)
            x1 .+= t
        end
        # Remove the solved part after all parts with equal seq number are dealt with.
        if d!=1<<length(bins[1])
            # b_ = b-A*x1
            apply!(A, b_, x1)
            b_ .= b .- b_
        end
    end
    x1
end

function timed_azs_solve(fplatform::BasisFunctions.Platform, i, f::Function; afirst=false, verbose=false, options...)
    t1 = @timed begin
        platform = fplatform.super_platform
        a = A(fplatform, i)
        zt = Zt(fplatform, i)
        s = sampler(fplatform, i)
        omega = grid(s)
        gamma = supergrid(omega)
        domain = FrameFun.domain(src(a))
        primal = FrameFun.basis(src(a))
        dual = BasisFunctions.wavelet_dual(primal)
        b = s*f
    end
    cutoff = default_cutoff(a)
    t2 = @timed bound = FrameFun.boundary_grid(gamma, domain)
    t3 = @timed coefficient_mask = BasisFunctions.coefficient_index_mask_of_overlapping_elements(dual, bound)
    t4 = @timed begin
        grid_mask = BasisFunctions.grid_index_mask_in_element_support(primal, gamma, coefficient_mask)
        grid_mask .= grid_mask .& FrameFun.mask(omega)
        DMZ = gamma[grid_mask]
        system_coefficient_mask = BasisFunctions.coefficient_index_mask_of_overlapping_elements(primal, DMZ)
        gr = restriction_operator(gridbasis(omega, coeftype(primal)), gridbasis(DMZ, coeftype(primal)))
        dr = restriction_operator(primal, system_coefficient_mask)
    end
    if afirst
        t5 = @timed m = matrix(gr*a*dr')
        t6 = @timed x2 = dr'*LAPACK.gelsy!(m,gr*b,cutoff)[1]
        t7 = @timed x1 = zt*(b-a*x2)
    else
        t7 = @timed x1 = zt*b
        t5 = @timed m = matrix(gr*a*dr')
        t6 = @timed x2 = dr'*LAPACK.gelsy!(m,gr*(b-a*x1),cutoff)[1]
    end
    t8 = @timed x1 .= x1 .+ x2
    info("error is $(norm(a*x1-b))")
    x1, [t1,t2,t3,t4,t5,t6,t7,t8]
end
