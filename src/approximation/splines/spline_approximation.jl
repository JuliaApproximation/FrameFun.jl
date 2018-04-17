
# Some of the util functions situated here can be located to more suitable locations.


# The main function
"""
The operators used for spline approximation using the AZ algorithm.

The domain, sizes and the degrees of the tensor spline frame is required.

"""
function spline_approximation_operators end

"""
Index of elements of `B` that overlap with `boundary`.
"""
function boundary_element_indices(B, boundary::MaskedGrid)
    s = Set{Int}()
    for x in boundary
        push!(s,overlapping_elements(B,x)...)
    end
    collect(s)
end


"""
A grid that contains the points of `omega_grid` that are not evaluated to zero by the elements that overlap with boundary_grid.
"""
function boundary_support_grid(B, boundary_grid::MaskedGrid, omega_grid::MaskedGrid)
    boundary_indices = boundary_element_indices(B,boundary_grid)
    s = Set{Int}()
    for i in boundary_indices
        push!(s,BasisFunctions.support_indices(B,supergrid(omega_grid),i)...)
    end
    a = collect(s)
    mask = zeros(Bool,size(omega_grid.mask))
    mask[a] = true
    mask .= mask .& omega_grid.mask
    MaskedGrid(supergrid(omega_grid),mask)
end

function spline_approximation_operators(domain, N1, N2, degr1, degr2)
    B1 = BSplineTranslatesBasis(N1,degr1);
    g1 = BasisFunctions.oversampled_grid(B1,2)
    B2 = BSplineTranslatesBasis(N2,degr2);
    g2 = BasisFunctions.oversampled_grid(B2,2)
    B = B1⊗B2;
    g = g1×g2;

    E1 = CirculantOperator(evaluation_matrix(B1[1],g1)[:])*IndexExtensionOperator(Span(B1),gridspace(g1),1:2:length(g1))
    E2 = CirculantOperator(evaluation_matrix(B2[1],g2)[:])*IndexExtensionOperator(Span(B2),gridspace(g2),1:2:length(g2))

    G1 = CirculantOperator(E1'E1*[1,zeros(length(g1)-1)...]);
    G2 = CirculantOperator(E2'E2*[1,zeros(length(g2)-1)...]);
    G = G1⊗G2
    E = E1⊗E2

    omega_grid = FrameFun.subgrid(g, domain);
    S = evaluation_operator(Span(B), omega_grid).operators[end]

    delta_omega_grid = boundary_grid(g, domain);
    SO2B = restriction_operator(omega_grid, delta_omega_grid);

    boundary_support = boundary_support_grid(B, delta_omega_grid, omega_grid);
    SO2S = restriction_operator(omega_grid, boundary_support);

    BEO = boundary_extension_operator(delta_omega_grid, B);

    DG = inv(G);
    Eo = S*E
    DEo = Eo*DG

    A = Eo;
    Z = DEo;

    I = IdentityOperator(dest(A));
    P = I-A*Z';

    Sampler = GridSamplingOperator(gridspace(omega_grid))

    Sampler, A, Z, P, SO2S, BEO, SO2B
end



function spline_approximation_operators(domain, N1, N2, N3, degr1, degr2, degr3)
    B1 = BSplineTranslatesBasis(N1,degr1);
    g1 = BasisFunctions.oversampled_grid(B1,2)
    B2 = BSplineTranslatesBasis(N2,degr2);
    g2 = BasisFunctions.oversampled_grid(B2,2)
    B3 = BSplineTranslatesBasis(N3,degr3);
    g3 = BasisFunctions.oversampled_grid(B3,2)
    B = B1⊗B2⊗B3 ;
    g = g1×g2×g3;

    E1 = CirculantOperator(evaluation_matrix(B1[1],g1)[:])*IndexExtensionOperator(Span(B1),gridspace(g1),1:2:length(g1))
    E2 = CirculantOperator(evaluation_matrix(B2[1],g2)[:])*IndexExtensionOperator(Span(B2),gridspace(g2),1:2:length(g2))
    E3 = CirculantOperator(evaluation_matrix(B3[1],g3)[:])*IndexExtensionOperator(Span(B3),gridspace(g3),1:2:length(g3))

    G1 = CirculantOperator(E1'E1*[1,zeros(length(g1)-1)...]);
    G2 = CirculantOperator(E2'E2*[1,zeros(length(g2)-1)...]);
    G3 = CirculantOperator(E3'E3*[1,zeros(length(g3)-1)...]);
    G = G1⊗G2⊗G3
    E = E1⊗E2⊗E3

    omega_grid = FrameFun.subgrid(g, domain);
    S = evaluation_operator(Span(B), omega_grid).operators[end]

    delta_omega_grid = boundary_grid(g, domain);
    SO2B = restriction_operator(omega_grid, delta_omega_grid);

    boundary_support = boundary_support_grid(B, delta_omega_grid, omega_grid);
    SO2S = restriction_operator(omega_grid, boundary_support);

    BEO = boundary_extension_operator(delta_omega_grid, B);

    DG = inv(G);
    Eo = S*E
    DEo = Eo*DG

    A = Eo;
    Z = DEo;

    I = IdentityOperator(dest(A));
    P = I-A*Z';

    Sampler = GridSamplingOperator(gridspace(omega_grid))

    Sampler, A, Z, P, SO2S, BEO, SO2B
end
