
using FrameFun.Platforms, FrameFun.FrameFunInterface, FrameFun.ExtensionFramePlatforms,
    FrameFun.FrameFunInterface, FrameFun.ExtensionFramePlatforms, FrameFun.ExtensionFrames,
    FrameFun.WeightedSumPlatforms, FrameFun.AugmentationPlatforms, FrameFun.FrameFunInterfaceExtension,
    BasisFunctions, FrameFun.ApproximationProblems,
    DomainSets, StaticArrays, Test, LinearAlgebra

# TODO add AugmentationPlatforms tests
d1 = 0.0..0.5
    plat1 = platform(extensionframe(Fourier(10),0.0..0.5)), 10
    d2 = (0.0..0.5)^2
    plat2 = platform(extensionframe(Fourier(10)^2,d2)), (10,10)
    d3 = disk(.45,SVector(.5,.5))
    plat3 = platform(extensionframe(Fourier(10)^2,d3)), (10,10)
    plat4 = platform(extensionframe(Fourier(10,-1,1)^2,d3)), (10,10)

    wplat1 = WeightedSumPlatform(platform(Fourier(10)), x->1., sqrt), (10)
    wplat2 = WeightedSumPlatform(platform(Fourier(10)^2), (x,y)->1., (x,y)->sqrt(x^2+y^2)), (10,10)


    ap1 = approximationproblem(Fourier(10),d1)
    ap2 = approximationproblem(Fourier(10)^2,d2)
    ap3 = approximationproblem(Fourier(10)^2,d3)
    ap4 = approximationproblem(Fourier(10,-1,1)^2,d3)
    wap1 = approximationproblem(wplat1...)
    wap2 = approximationproblem(wplat2...)

    dict1 = dictionary(ap1)
    dict2 = dictionary(ap2)
    dict3 = dictionary(ap3)
    dict4 = dictionary(ap4)
    wdict1 = dictionary(wap1)
    wdict2 = dictionary(wap2)

@testset "dictionary" begin
    for plat in (plat1,plat2,plat3,plat4)
        @test dictionary(plat...) isa ExtensionFrameSuper
        @test dualdictionary(plat...,measure(plat[1])) isa ExtensionFrameSuper
    end

    for ap in (ap1,ap2,ap3,ap4)
        @test dictionary(ap) isa ExtensionFrameSuper
    end

    for (plat,ap) in zip((wplat1,wplat2),(wap1,wap2))
        @test dictionary(plat...) isa MultiDict
        @test dictionary(ap) isa MultiDict
        @test dualdictionary(plat...,measure(plat[1])) isa MultiDict
    end
end


@testset "samplingparameter" begin
    @test samplingparameter(ap1) == 38
    @test samplingparameter(ap2) == (27,27)
    @test samplingparameter(ap3) == (18,18)
    @test samplingparameter(ap4) == (36,36)

    @test samplingparameter(wap1) == (40)
    @test samplingparameter(wap2) == (21,21)


    @test samplingparameter(plat1...)  == samplingparameter(ap1)
    @test samplingparameter(plat2...)   == samplingparameter(ap2)
    @test samplingparameter(plat3...) == samplingparameter(ap3)
    @test samplingparameter(plat4...) == samplingparameter(ap4)
    @test samplingparameter(wplat1...) == samplingparameter(wap1)
    @test samplingparameter(wplat2...) == samplingparameter(wap2)
end

@testset "oversampling_grid" begin

    @test oversampling_grid(ap1) isa IndexSubGrid
    @test oversampling_grid(ap2) isa TensorSubGrid
    @test oversampling_grid(ap3) isa GridArrays.MaskedGrid
    @test oversampling_grid(ap4) isa GridArrays.MaskedGrid

    @test oversampling_grid(wap1) == FourierGrid(40)
    @test oversampling_grid(wap2) == FourierGrid(21)^2

    @test oversampling_grid(plat1...) isa IndexSubGrid
    @test oversampling_grid(plat2...) isa TensorSubGrid
    @test oversampling_grid(plat3...) isa GridArrays.MaskedGrid
    @test oversampling_grid(plat4...) isa GridArrays.MaskedGrid
    @test oversampling_grid(wplat1...) == FourierGrid(40)
    @test oversampling_grid(wplat2...) == FourierGrid(21)^2
end

@testset "samplingoperator" begin
    f1 = exp
    f2 = (x,y) -> exp(x*y)
    op1 = samplingoperator(ap1)
    op2 = samplingoperator(ap2)
    op3 = samplingoperator(ap3)
    op4 = samplingoperator(ap4)

    op5 = samplingoperator(wap1)
    op6 = samplingoperator(wap2)
    # Generic operators so hard to check their structure

    @test samplingoperator(ap1)*f1≈samplingoperator(plat1...)*f1
    @test samplingoperator(ap2)*f2≈samplingoperator(plat2...)*f2
    @test samplingoperator(ap3)*f2≈samplingoperator(plat3...)*f2
    @test samplingoperator(ap4)*f2≈samplingoperator(plat4...)*f2

    @test samplingoperator(wap1)*f1≈samplingoperator(wplat1...)*f1
    @test samplingoperator(wap2)*f2≈samplingoperator(wplat2...)*f2
end


@testset "sampling_grid" begin
    for ap in (ap1,ap2,ap3,ap4,wap1,wap2)
        @test sampling_grid(ap)≈oversampling_grid(ap)
    end
    for plat in (plat1,plat2,plat3,plat4,wplat1,wplat2)
        @test sampling_grid(plat...)≈oversampling_grid(plat...)
    end
    for dict in (dict1,dict2,dict3,dict4)
        @test sampling_grid(dict)≈oversampling_grid(dict)
    end
end

@testset "discretemeasure" begin
    for ap in (ap1,ap2,ap3,ap4,wap1,wap2)
        @test discretemeasure(ap)≈discretemeasure(sampling_grid(ap))
    end
    for plat in (plat1,plat2,plat3,plat4,wplat1,wplat2)
        @test discretemeasure(plat...)≈discretemeasure(sampling_grid(plat...))
    end
    for dict in (dict1,dict2,dict3,dict4)
        @test discretemeasure(dict)≈discretemeasure(sampling_grid(dict))
    end
end


@testset "measure" begin
    @test measure(ap1)==measure(dict1)
    @test measure(ap2)==measure(dict2)
    @test measure(ap3)==measure(dict3)
    @test measure(ap4)==measure(dict4)
    @test measure(wap1)==measure(wdict1)
    @test measure(wap2)==measure(wdict2)
    @test measure(plat1...)==measure(dict1)
    @test measure(plat2...)==measure(dict2)
    @test measure(plat3...)==measure(dict3)
    @test measure(plat4...)==measure(dict4)
    @test measure(wplat1...)==measure(wdict1)
    @test measure(wplat2...) ==measure(wdict2)
end


@testset "azdual_dict" begin
    @test azdual_dict(ap1) isa ExtensionFrame
    @test azdual_dict(ap2) isa ExtensionFrameTensor
    @test azdual_dict(ap3) isa ExtensionFrame
    @test azdual_dict(ap4) isa ExtensionFrame
    @test azdual_dict(wap1) isa MultiDict
    @test azdual_dict(wap2) isa MultiDict

    @test azdual_dict(plat1...) isa ExtensionFrame
    @test azdual_dict(plat2...) isa ExtensionFrameTensor
    @test azdual_dict(plat3...) isa ExtensionFrame
    @test azdual_dict(plat4...) isa ExtensionFrame
    @test azdual_dict(wplat1...) isa MultiDict
    @test azdual_dict(wplat2...) isa MultiDict

    @test azdual_dict(dict1) isa ExtensionFrame
    @test azdual_dict(dict2) isa ExtensionFrameTensor
    @test azdual_dict(dict3) isa ExtensionFrame
    @test azdual_dict(dict4) isa ExtensionFrame
end

@testset "discretization" begin
    op = discretization(ap1)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict1, r, sampling_grid(ap1))
    op = discretization(ap2)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict2, r, sampling_grid(ap2))
    op = discretization(ap3)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict3, r, sampling_grid(ap3))
    op = discretization(ap4)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict4, r, sampling_grid(ap4))
    op = discretization(wap1)
    r = rand(src(op))
    @test op*r≈eval_expansion(wdict1, r, sampling_grid(wap1))
    op = discretization(wap2)
    r = rand(src(op))
    @test op*r≈eval_expansion(wdict2, r, sampling_grid(wap2))

    op = discretization(plat1...)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict1, r, sampling_grid(ap1))
    op = discretization(plat2...)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict2, r, sampling_grid(ap2))
    op = discretization(plat3...)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict3, r, sampling_grid(ap3))
    op = discretization(plat4...)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict4, r, sampling_grid(ap4))
    op = discretization(wplat1...)
    r = rand(src(op))
    @test op*r≈eval_expansion(wdict1, r, sampling_grid(wap1))
    op = discretization(wplat2...)
    r = rand(src(op))
    @test op*r≈eval_expansion(wdict2, r, sampling_grid(wap2))

    op = discretization(dict1)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict1, r, sampling_grid(ap1))
    op = discretization(dict2)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict2, r, sampling_grid(ap2))
    op = discretization(dict3)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict3, r, sampling_grid(ap3))
    op = discretization(dict4)
    r = rand(src(op))
    @test op*r≈eval_expansion(dict4, r, sampling_grid(ap4))

end

@testset "dualdiscretization" begin
    op = dualdiscretization(ap1)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap1), r, sampling_grid(ap1))
    op = dualdiscretization(ap2)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap2), r, sampling_grid(ap2))
    op = dualdiscretization(ap3)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap3), r, sampling_grid(ap3))
    op = dualdiscretization(ap4)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap4), r, sampling_grid(ap4))
    op = dualdiscretization(wap1)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(wap1), r, sampling_grid(wap1))
    op = dualdiscretization(wap2)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(wap2), r, sampling_grid(wap2))

    op = dualdiscretization(plat1...)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap1), r, sampling_grid(ap1))
    op = dualdiscretization(plat2...)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap2), r, sampling_grid(ap2))
    op = dualdiscretization(plat3...)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap3), r, sampling_grid(ap3))
    op = dualdiscretization(plat4...)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap4), r, sampling_grid(ap4))
    op = dualdiscretization(wplat1...)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(wap1), r, sampling_grid(wap1))
    op = dualdiscretization(wplat2...)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(wap2), r, sampling_grid(wap2))

    op = dualdiscretization(dict1)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap1), r, sampling_grid(ap1))
    op = dualdiscretization(dict2)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap2), r, sampling_grid(ap2))
    op = dualdiscretization(dict3)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap3), r, sampling_grid(ap3))
    op = dualdiscretization(dict4)
    r = rand(src(op))
    @test op*r≈eval_expansion(azdual_dict(ap4), r, sampling_grid(ap4))
end


@testset "solvers" begin
    for dict in (dict1,)
        op1 = solver(dict;solverstyle=DirectStyle())
        r = rand(src(op1))
        op2 = solver(dict;solverstyle=IterativeStyle())
        @test norm(op1*r-op2*r) < 1e-4
        op3 = solver(dict;solverstyle=AZStyle())
        @test norm(op1*r-op3*r) < 1e-4
        op4 = solver(dict;solverstyle=AZSmoothStyle())
        @test norm(op1*r-op4*r) < 1e-4
        op4 = solver(dict;solverstyle=TridiagonalProlateStyle())
    end

    for dict in (dict1,dict3,dict4)
        op1 = solver(dict;solverstyle=DirectStyle())
        r = rand(src(op1))
        op2 = solver(dict;solverstyle=AZStyle())
        @test  norm(op1*r-op2*r) < 1e-1
    end

    P = FourierExtensionPlatform(0.0..0.5)
    f = exp
    A,b = discretization(f, P, 100)
    M = plungeoperator(P,100)*A
    S = rSVD_solver(M;threshold=1e-4)
    x1 = S*plungeoperator(P,100)*b
    x2 = AZ_Zt(P,100)*(b-A*x1)
    @test norm(A*(x1+x2)-b) <  1e-3
end

@testset "AZ_A, AZ_Z, AZ_Zt, plungeoperator, firstAZstepoperator, plungerank" begin

    for (ap, plat) in ((ap1,plat1),(ap2,plat2),(ap3,plat3),(ap4,plat4),(wap1,wplat1),(wap2,wplat2))
        @test AZ_A(ap) ≈ AZ_A(plat...)
        @test AZ_Zt(ap) ≈ AZ_Zt(plat...)
        @test plungeoperator(ap) ≈ plungeoperator(plat...)
        @test firstAZstepoperator(ap) ≈ firstAZstepoperator(plat...)
        @test abs(plungerank(ap) - plungerank(plat...)) < 3
    end

    for (ap, dict) in ((ap1,dict1),(ap2,dict2),(ap3,dict3),(ap4,dict4),)
        @test AZ_A(ap) ≈ AZ_A(dict)
        @test AZ_Zt(ap) ≈ AZ_Zt(dict)
        @test plungeoperator(ap) ≈ plungeoperator(dict)
        @test firstAZstepoperator(ap) ≈ firstAZstepoperator(dict)
        @test abs(plungerank(ap) - plungerank(dict)) < 3
    end

    for ap in (ap1,ap2,ap3,ap4,wap1,wap2)
        @test norm(plungeoperator(ap) - (I -AZ_A(ap)*AZ_Zt(ap))) < 1e-12
        @test norm(firstAZstepoperator(ap) - (AZ_A(ap) -AZ_A(ap)*AZ_Zt(ap)*AZ_A(ap) )) < 1e-12
    end

    @test plungerank(plat1[1],(100);threshold=1e-6) <= 30
    @test plungerank(plat2[1],(30,30);threshold=1e-6) <= 700
    @test plungerank(plat3[1],(30,30);threshold=1e-6) <= 600
    @test plungerank(plat4[1],(30,30);threshold=1e-6) <= 600
    @test plungerank(wplat1[1],100;threshold=1e-6) <= 30
    @test plungerank(wplat2[1],(15,15);threshold=1e-6)  <= 400

end

@testset "approximate" begin
    P = FourierExtensionPlatform(0.0..0.5)
    f = exp
    A,b = discretization(f, P, 100)
    M = plungeoperator(P,100)*A
    S = rSVD_solver(M;threshold=1e-4)
    x1 = S*plungeoperator(P,100)*b
    x2 = AZ_Zt(P,100)*(b-A*x1)
    @test norm(A*(x1+x2)-b) <  1e-3

    F = approximate(f, P, 100;REG=rSVD_solver)[1]
    x = .2; @test abs(F(x)-f(x)) < 1e-3

    P = FourierExtensionPlatform(0.0..0.5)
    f = exp
    A,b = discretization(f, P, 100)
    M = plungeoperator(P,100)*A
    S = pQR_solver(M;threshold=1e-4)
    x1 = S*plungeoperator(P,100)*b
    x2 = AZ_Zt(P,100)*(b-A*x1)
    @test norm(A*(x1+x2)-b) < 1e-3

    F = approximate(f, P, 100;REG=pQR_solver)[1]
    x = .2; @test abs(F(x)-f(x)) < 1e-3

    P = FourierExtensionPlatform(0.0..0.5)
    f = exp
    A,b = discretization(f, P, 100)
    M = plungeoperator(P,100)*A
    S = pSVD_solver(M;threshold=1e-4)
    x1 = S*plungeoperator(P,100)*b
    x2 = AZ_Zt(P,100)*(b-A*x1)
    @test norm(A*(x1+x2)-b) < 1e-3

    F = approximate(f, P, 100;REG=pSVD_solver)[1]
    x = .2; @test abs(F(x)-f(x)) < 1e-3
end
