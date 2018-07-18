module test_suite_applications
using BasisFunctions
using FrameFun
using FrameFun: boundary_support_grid, relative_indices, restriction_operator, azselection_restriction_operators
using WaveletsCopy: cdf24, db3
using Base.Test
using StaticArrays
using Domains
# Select the frame functions that overlap with the boundary


function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end
delimit("Spline approximation")
@testset "Util" begin
    # center = @SVector [.5,.5]
    # domain  = disk(.3,center)
    # N1 = 40; N2 = 40;
    # degr1 = 2; degr2 = 3;
    degree = [2,3]
    init = [40,40]
    i = 1
    oversampling = 2
    center = @SVector [.5,.5]
    domain  = disk(.3,center)
    T = Float64
    platform=bspline_platform(T, init, degree, oversampling)
    fplatform = extension_frame_platform(platform, domain)

    B = superdict(primal(fplatform, i))
    S, R = azselection_restriction_operators(fplatform, i)
    @test size(S)==(476, 1600)
    @test size(R)==(988 , 1808)

    omega_grid = BasisFunctions.grid(sampler(fplatform, 1))
    g = supergrid(omega_grid)

    g1 = boundary_grid(g, domain);
    g2 = boundary_support_grid(B,g1,omega_grid)

    @test relative_indices(g1,omega_grid)==[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 31, 32, 33, 34, 35, 36, 37, 56, 57, 58, 59, 60, 83, 84, 85, 86, 87, 112, 113, 114, 115, 116, 145, 146, 147, 148, 179, 180, 181, 182, 215, 216, 217, 218, 253, 254, 255, 292, 293, 294, 331, 332, 333, 334, 373, 374, 375, 416, 417, 418, 459, 460, 461, 504, 505, 506, 549, 550, 551, 596, 597, 642, 643, 688, 689, 690, 735, 736, 737, 784, 785, 832, 833, 880, 881, 928, 929, 976, 977, 1024, 1025, 1072, 1073, 1074, 1119, 1120, 1121, 1166, 1167, 1212, 1213, 1258, 1259, 1260, 1303, 1304, 1305, 1348, 1349, 1350, 1391, 1392, 1393, 1434, 1435, 1436, 1475, 1476, 1477, 1478, 1515, 1516, 1517, 1554, 1555, 1556, 1591, 1592, 1593, 1594, 1627, 1628, 1629, 1630, 1661, 1662, 1663, 1664, 1693, 1694, 1695, 1696, 1697, 1722, 1723, 1724, 1725, 1726, 1749, 1750, 1751, 1752, 1753, 1772, 1773, 1774, 1775, 1776, 1777, 1778, 1791, 1792, 1793, 1794, 1795, 1796, 1797, 1798, 1799, 1800, 1801, 1802, 1803, 1804, 1805, 1806, 1807, 1808]

    @test size(restriction_operator(gridbasis(omega_grid),gridbasis(g1)))==(186, 1808)
end

@testset "Spline approximation" begin
    init = [3,3]
    degree = [1,3]
    oversampling = 2
    center = @SVector [.5,.5]
    domain  = disk(.3,center)
    i = 2
    # correct caclulation of the A and Z matrix
    for T in [Float64, BigFloat]
        center = @SVector [T(.5),T(.5)]
        domain  = disk(T(.3),center)
        platform=BasisFunctions.bspline_platform(T, init, degree, oversampling)
        fplatform = extension_frame_platform(platform, domain)

        P = primal(fplatform,i)
        D = dual(fplatform,i)
        S = sampler(fplatform, i)

        EP = BasisFunctions.A(fplatform, i)
        ED = BasisFunctions.Zt(fplatform, i)'

        p = fplatform.parameter_sequence[i]
        Pt = tensorproduct([BSplineTranslatesBasis(pi, di, T) for (pi, di) in zip(p, degree) ])
        g = BasisFunctions.oversampled_grid(Pt, oversampling)
        Ft = extensionframe(Pt, domain)
        gO = FrameFun.subgrid(g, domain)
        EPt = evaluation_operator(Pt, gO)
        EDt = EPt*DiscreteDualGram(Pt,oversampling=oversampling)
        e = map(T,rand(size(Pt)...))
        @test EPt*e≈EP*e
        @test EDt*e≈ED*e*length(Pt)
    end

    # Correct implementation of the az and azs Algorithm
    # for irregular domain
    f2d = (x,y) -> x*(y-1)^2
    center = @SVector [.5,.5]
    domain2d = disk(.3,center)
    epsilon=1e-10
    degree = [1,2]
    init = [20,20]
    oversampling = 2
    platform = bspline_platform(Float64, init, degree, oversampling)
    fplatform = extension_frame_platform(platform, domain2d)
    i = 1
    s = sampler(fplatform, i)
    a = BasisFunctions.A(fplatform, i)
    zt = BasisFunctions.Zt(fplatform, i)
    p = FrameFun.plunge_operator(a,zt)
    rd,sb = azselection_restriction_operators(fplatform, i)
    b = s*f2d
    r = FrameFun.estimate_plunge_rank(a)
    @test r==102
    AZ = AZSolver(a,zt,R=r,cutoff=epsilon)
    AZS = AZSSolver(a,zt,rd', sb,cutoff=epsilon)
    x = zeros(src(a))

    apply!(AZ, x, b)
    @test norm(a*x-b)+1≈ 1
    apply!(AZS, x, b)
    @test norm(a*x-b)+1≈ 1
    x = AZS*b
    @test norm(a*x-b)+1≈ 1
    x = FrameFun.az_solve(b, a, zt, cutoff=epsilon, R=r)
    @test norm(a*x-b)+1≈ 1
    x = FrameFun.az_solve(fplatform, i, f2d; cutoff=epsilon)
    @test norm(a*x-b)+1≈ 1
    x = FrameFun.azs_solve(fplatform, i, f2d; cutoff=epsilon)
    @test norm(a*x-b)+1≈ 1

    # For tensor domain
    i = 2
    platform  = bspline_platform(Float64, [4,4], [1,1], 2)
    fun = (x,y)->1+x+y
    mid = .25
    dom = interval(0,.5)^2
    fplatform = extension_frame_platform(platform, dom)
    p = primal(platform, i);
    s = sampler(platform, i);
    BR,DMZ_R = azselection_restriction_operators(fplatform, i)
    A = BasisFunctions.A(fplatform, i);
    Zt = BasisFunctions.Zt(fplatform, i);
    r = FrameFun.estimate_plunge_rank(A)
    @test r==12
    S = sampler(fplatform, i)
    boundary = FrameFun.boundary_grid(BasisFunctions.grid(s), dom)
    dx = BasisFunctions.stepsize(elements(BasisFunctions.grid(s))[1])
    split_domain = interval(mid-dx/2,mid+dx/2)×interval(0,1)
    split_grid = FrameFun.subgrid(boundary, split_domain);
    DMZsplit = FrameFun.boundary_support_grid(p, split_grid, boundary)
    BRsplit, DMZ_Rsplit = FrameFun._spline_util_restriction_operators(p, BasisFunctions.grid(S), DMZsplit)
    left_over = boundary-DMZsplit
    boundary1, boundary2 = split(left_over)
    BR1, DMZ_R1 = FrameFun._spline_util_restriction_operators(p, BasisFunctions.grid(S), boundary1)
    BR2, DMZ_R2 = FrameFun._spline_util_restriction_operators(p, BasisFunctions.grid(S), boundary2)


    b = S*fun
    AZSS = FrameFun.AZSSolver(A, Zt, BR', DMZ_R)
    x = AZSS*b
    @test 1+norm(A*x-b)≈1
    AZS = FrameFun.AZSolver(A, Zt)
    x = AZS*b
    @test 1+norm(A*x-b)≈1

    x = FrameFun.az_solve(b, A, Zt)
    @test 1+norm(A*x-b)≈1
    x = FrameFun.azs_solve(b, A, Zt, BR', DMZ_R)
    @test 1+norm(A*x-b)≈1

    x = FrameFun.az_solve(fplatform, i, fun)
    @test 1+norm(A*x-b)≈1
    x = FrameFun.azs_solve(fplatform, i, fun)
    @test 1+norm(A*x-b)≈1

    op = FrameFun.AZSSolver(fplatform, i)
    x = op*b;@test 1+norm(A*x-b)≈1
    op = FrameFun.AZSolver(fplatform, i)
    x = op*b;@test 1+norm(A*x-b)≈1

    # 1D
    f1d = identity
    domain1d = interval(0,.5)
    epsilon=1e-10
    degree = 1
    init = 20
    oversampling = 2
    platform = bspline_platform(Float64, init, degree, oversampling)
    fplatform = extension_frame_platform(platform, domain1d)
    i = 4
    s = sampler(fplatform, i)
    a = BasisFunctions.A(fplatform, i)
    zt = BasisFunctions.Zt(fplatform, i)
    p = FrameFun.plunge_operator(a,zt)
    rd,sb = azselection_restriction_operators(fplatform, i)
    b = s*f1d
    r = FrameFun.estimate_plunge_rank(a)
    @test r==2
    AZ = AZSolver(a,zt,R=r,cutoff=epsilon)
    AZS = AZSSolver(a,zt,rd', sb,cutoff=epsilon)
    x = zeros(src(a))

    apply!(AZ, x, b)
    @test norm(a*x-b)+1≈ 1
    apply!(AZS, x, b)
    @test norm(a*x-b)+1≈ 1
    x = AZS*b
    @test norm(a*x-b)+1≈ 1
    x = FrameFun.az_solve(b, a, zt, cutoff=epsilon, R=r)
    @test norm(a*x-b)+1≈ 1
    x = FrameFun.az_solve(fplatform, i, f1d; cutoff=epsilon)
    @test norm(a*x-b)+1≈ 1
    x = FrameFun.azs_solve(fplatform, i, f1d; cutoff=epsilon)
    @test norm(a*x-b)+1≈ 1
    # warn("az tree not in 1D")
    # x = FrameFun.az_tree_solve(fplatform, i, f1d; cutoff=epsilon)
    # @test norm(a*x-b)+1≈ 1
    op = FrameFun.AZSSolver(fplatform, i)
    x = op*b;@test 1+norm(a*x-b)≈1
    op = FrameFun.AZSolver(fplatform, i)
    x = op*b;@test 1+norm(a*x-b)≈1


end

@testset "Scaling approximation" begin
    # For tensor domain
    i = 2
    os = 1; fun = (x,y)->1.;
    platform = scaling_platform([4,5], [cdf24,cdf24], 1<<os)
    dom = interval(0,.5)^2
    fplatform = extension_frame_platform(platform, dom)
    p = primal(platform, i);
    s = sampler(platform, i);
    A = BasisFunctions.A(fplatform, i);
    Zt = BasisFunctions.Zt(fplatform, i);
    S = sampler(fplatform, i)
    BR,DMZ_R = azselection_restriction_operators(fplatform, i)

    b = S*fun
    AZSS = FrameFun.AZSSolver(A, Zt, BR', DMZ_R;afirst=false)
    x = AZSS*b
    @test 1+norm(A*x-b)≈1

    op = FrameFun.AZSSolver(fplatform, i;afirst=false)
    x = op*b;@test 1+norm(A*x-b)≈1

    x = FrameFun.azs_solve(b, A, Zt, BR', DMZ_R;afirst=false)
    @test 1+norm(A*x-b)≈1

    x = FrameFun.azs_solve(fplatform, i, fun;afirst=false)
    @test 1+norm(A*x-b)≈1

    # # 1D
    f1d = identity
    domain1d = interval(0,.5)
    epsilon=1e-10
    platform = scaling_platform(4, db3, 1<<os)
    fplatform = extension_frame_platform(platform, domain1d)
    i = 2
    s = sampler(fplatform, i)
    a = BasisFunctions.A(fplatform, i)
    zt = BasisFunctions.Zt(fplatform, i)
    rd,sb = azselection_restriction_operators(fplatform, i)
    b = s*f1d
    AZS = AZSSolver(a,zt,rd', sb,cutoff=epsilon,afirst=false)
    x = zeros(src(a))

    apply!(AZS, x, b)
    @test norm(a*x-b)+1≈ 1
    x = AZS*b
    @test norm(a*x-b)+1≈ 1
    x = FrameFun.azs_solve(fplatform, i, f1d; cutoff=epsilon,afirst=false)
    @test norm(a*x-b)+1≈ 1
    op = FrameFun.AZSSolver(fplatform, i, afirst=false)
    x = op*b;@test 1+norm(a*x-b)≈1

end
end
