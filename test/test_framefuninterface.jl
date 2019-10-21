using FrameFun.FrameFunInterface, FrameFun.Platforms, FrameFun.ApproximationProblems,
    Test, LinearAlgebra, BasisFunctions, FrameFun.ParameterPaths, FrameFun.WeightedSumPlatforms,
    FrameFun.ExtensionFramePlatforms


ap1 = approximationproblem(platform(Fourier(10)),10)
ap2 = approximationproblem(platform(Fourier(10)^2),(10,10))
dict1 = Fourier(10)
dict2 = Fourier(10)^2
plat1 = (platform(Fourier(10)),10)
plat2 = (platform(Fourier(10)^2),(10,10))

aps = (ap1,ap2,)
alls = (dict1, dict2, aps...,)

@testset "samplingparameter" begin
    # adaptive approximations has no samplingparameter
    for (a,sol) in zip(alls ,(10,(10,10),10,(10,10)))
        @test samplingparameter(a)==sol
    end
    @test samplingparameter(plat1...) == 10
    @test samplingparameter(plat2...) == (10,10)
end

@testset "interpolation_grid" begin
    @test interpolation_grid(plat1...) == interpolation_grid(dict1)
    @test interpolation_grid(plat2...) == interpolation_grid(dict2)

    @test interpolation_grid(ap1) == interpolation_grid(dict1)
    @test interpolation_grid(ap2) == interpolation_grid(dict2)
end

@testset "oversampling_grid" begin
    @test oversampling_grid(dict1) == FourierGrid(10)
    @test oversampling_grid(dict1, 11) == FourierGrid(11)
    @test oversampling_grid(dict2) == FourierGrid(10)^2
    @test oversampling_grid(dict2, (11,11)) == FourierGrid(11)^2
    @test oversampling_grid(plat1...) == FourierGrid(10)
    @test oversampling_grid(ap1) == FourierGrid(10)
    @test oversampling_grid(ap2) == FourierGrid(10)^2
end

@testset "samplingoperator" begin
    op1 = samplingoperator(dict1)
    op2 = samplingoperator(dict2)
    op3 = samplingoperator(plat1...)
    op4 = samplingoperator(plat2...)
    op5 = samplingoperator(ap1)
    op6 = samplingoperator(ap2)

    op1a, op1b = elements(op1)
    @test grid(op1a) == FourierGrid(10)
    @test op1b ≈ IdentityOperator(dict1)

    op2a, op2b = elements(op3)
    @test grid(op2a) == FourierGrid(10)
    @test op2b ≈ IdentityOperator(dict1)

    op3a, op3b = elements(op5)
    @test grid(op3a) == FourierGrid(10)
    @test op3b ≈ IdentityOperator(dict1)
    f2 = (x,y)->exp(x*y)
    @test tensorproduct(op1,op1)*f2≈op2*f2
    @test tensorproduct(op3,op3)*f2≈op4*f2
    @test tensorproduct(op5,op5)*f2≈op6*f2
end


@testset "sampling_grid" begin
    @test sampling_grid(dict1) == FourierGrid(10)
    @test sampling_grid(dict2) == FourierGrid(10)^2
    @test sampling_grid(plat1...) == FourierGrid(10)
    @test sampling_grid(plat2...) == FourierGrid(10)^2
    @test sampling_grid(ap1) == FourierGrid(10)
    @test sampling_grid(ap2) == FourierGrid(10)^2
end


@testset "platform_grid" begin
    @test platform_grid(dict1;samplingstyle=GridStyle(),grid=FourierGrid(1)) == FourierGrid(1)
    @test platform_grid(dict2;samplingstyle=ProductSamplingStyle(GridStyle(),GridStyle()),grid=FourierGrid(1)) == FourierGrid(1)^2
    @test platform_grid(dict2;samplingstyle=GridStyle(),grid=FourierGrid(1)) == FourierGrid(1)
    @test platform_grid(plat1...,grid=FourierGrid(1)) == FourierGrid(1)
    @test platform_grid(plat2...,grid=FourierGrid(1)) == FourierGrid(1)^2
    @test platform_grid(ap1;samplingstyle=GridStyle(),grid=FourierGrid(1)) == FourierGrid(1)
    @test platform_grid(ap2;samplingstyle=GridStyle(),grid=FourierGrid(1)) == FourierGrid(1)
    @test sampling_grid(dict2;samplingstyle=GridStyle(),grid=FourierGrid(1))==FourierGrid(1)
    @test sampling_grid(dict2;samplingstyle=ProductSamplingStyle(GridStyle(),GridStyle()),grid=FourierGrid(1)) == FourierGrid(1)^2
end

@testset "discretemeasure" begin
    @test discretemeasure(dict1) == discretemeasure(sampling_grid(dict1))
    @test discretemeasure(dict2)== discretemeasure(sampling_grid(dict2))
    @test discretemeasure(plat1...)== discretemeasure(sampling_grid(plat1...))
    @test discretemeasure(plat2...)== discretemeasure(sampling_grid(plat2...))
    @test discretemeasure(ap1)== discretemeasure(sampling_grid(ap1))
    @test discretemeasure(ap2)== discretemeasure(sampling_grid(ap2))
end


@testset "measure" begin
    @test measure(dict1) == measure(dict1)
    @test measure(dict2)== measure(dict2)
    @test measure(plat1...)== measure(dict1)
    @test measure(plat2...)== measure(dict2)
    @test measure(ap1)== measure(dict1)
    @test measure(ap2)== measure(dict2)
end


@testset "azdual_dict" begin
    @test operator(azdual_dict(dict1))≈operator(azdual_dict(ap1))≈operator(azdual_dict(plat1...))
    op = azdual_dict(dict2)
    @test op isa TensorProductDict
    @test ≈(operator.(elements(op))...)
    op = azdual_dict(plat2...)
    @test op isa TensorProductDict
    @test ≈(operator.(elements(op))...)
    op = azdual_dict(ap2)
    @test op isa TensorProductDict
    @test ≈(operator.(elements(op))...)
    @test operator(element(azdual_dict(dict2),1))≈operator(azdual_dict(dict1))

    @test operator(azdual_dict(dict1;samplingstyle=GramStyle()))≈
        operator(element(azdual_dict(ap2;samplingstyle=ProductSamplingStyle(GramStyle(),GramStyle())),1))
    @test operator(azdual_dict(dict1;samplingstyle=GramStyle()))≈
        operator(element(azdual_dict(plat2...;samplingstyle=ProductSamplingStyle(GramStyle(),GramStyle())),1))
    @test operator(azdual_dict(dict1;samplingstyle=GramStyle()))≈
        operator(element(azdual_dict(dict2;samplingstyle=ProductSamplingStyle(GramStyle(),GramStyle())),1))
end

@testset "discretization" begin
    op1 = discretization(ap1)
    op2 = discretization(ap2)
    op3 = discretization(plat1...)
    op4 = discretization(plat2...)
    op5 = discretization(dict1)
    op6 = discretization(dict2)

    @test op2 isa TensorProductOperator
    @test op4 isa TensorProductOperator
    @test op6 isa TensorProductOperator

    @test op1 ≈ op3 ≈ op5
    @test op2 ≈ op4 ≈ op6

    f1 = exp
    f2 = (x,y)->exp(x*y)
    op1, b1 = discretization(f1,ap1)
    op2, b2 = discretization(f2,ap2)
    op3, b3 = discretization(f1,plat1...)
    op4, b4 = discretization(f2,plat2...)
    op5, b5 = discretization(f1,dict1)
    op6, b6 = discretization(f2,dict2)

    @test op2 isa TensorProductOperator
    @test op4 isa TensorProductOperator
    @test op6 isa TensorProductOperator

    @test op1 ≈ op3 ≈ op5
    @test op2 ≈ op4 ≈ op6

    @test b1≈b3≈b5
    @test b2≈b4≈b6
end

@testset "dualdiscretization" begin

    op1 = dualdiscretization(ap1)
    op2 = dualdiscretization(ap2)
    op3 = dualdiscretization(plat1...)
    op4 = dualdiscretization(plat2...)
    op5 = dualdiscretization(dict1)
    op6 = dualdiscretization(dict2)

    @test op2 isa TensorProductOperator
    @test op4 isa TensorProductOperator
    @test op6 isa TensorProductOperator

    @test op1 ≈ op3 ≈ op5
    @test op2 ≈ op4 ≈ op6

end

@testset "solver" begin
    op1 = solver(dict1)
    op2 = solver(dict2)
    op3 = solver(plat1...)
    op4 = solver(plat2...)
    op5 = solver(ap1)
    op6 = solver(ap2)

    @test op2 isa TensorProductOperator
    @test op4 isa TensorProductOperator
    @test op6 isa TensorProductOperator

    @test op1 ≈ op3 ≈ op5
    @test op2 ≈ op4 ≈ op6

    @test element(op2,1) ≈ op1 ≈ element(op2,2)


    for STYLE in (InverseStyle, DirectStyle, DualStyle, IterativeStyle,AZStyle,AZSmoothStyle)

        op1 = solver(dict1;solverstyle=STYLE())
        op2 = solver(dict2;solverstyle=ProductSolverStyle(STYLE(),STYLE()))
        op3 = solver(plat1...;solverstyle=STYLE())
        op4 = solver(plat2...;solverstyle=ProductSolverStyle(STYLE(),STYLE()))
        op5 = solver(ap1;solverstyle=STYLE())
        op6 = solver(ap2;solverstyle=ProductSolverStyle(STYLE(),STYLE()))



        @test op2 isa TensorProductOperator
        @test op4 isa TensorProductOperator
        @test op6 isa TensorProductOperator

        @test op1 ≈ op3 ≈ op5
        @test op2 ≈ op4 ≈ op6

        @test element(op2,1) ≈ op1 ≈ element(op2,2)
    end
end

@testset "AZ_A, AZ_Z, AZ_Zt, plungeoperator, firstAZstepoperator, plungerank" begin
    A1 = AZ_A(dict1)
    A2 = AZ_A(dict2)
    @test tensorproduct(A1,A1) ≈ A2

    Z1 = AZ_Z(dict1)
    Z2 = AZ_Z(dict2)
    @test tensorproduct(Z1,Z1) ≈ Z2

    Zt1 = AZ_Zt(dict1)
    Zt2 = AZ_Zt(dict2)
    @test tensorproduct(Zt1,Zt1) ≈ Zt2

    @test Zt2 ≈ Z2'

    P1 = plungeoperator(dict1)
    P2 = plungeoperator(dict2)
    @test 1+norm(tensorproduct(P1,P1) -P2)≈1
    @test I-A2*Zt2 ≈ P2

    M1 = firstAZstepoperator(dict1)
    M2 = firstAZstepoperator(dict2)
    @test 1+norm(tensorproduct(M1,M1) - M2)≈1
    @test 1+norm(A2*P2 - M2)≈1

    @test plungerank(dict1) <= 1
    @test plungerank(dict2) <= 1
end

@testset begin "ParameterPath"
    p = IncrementalCartesianParameterPath{2}()
    P = ExtensionFramePlatform(WeightedSumPlatform(platform(ChebyshevT(10,-1,1)^2), (x,y)->1.,
                (x,y)->sqrt(x^2+y^2)),.9*UnitDisk())
    path = HierarchyPath(p,ProductPath(p,p))
    paramP = parametrizedplatform(P, path)

    paramPparam=45
    Pparam = path[paramPparam]

    @test SamplingStyle(P)==SamplingStyle(paramP)
    @test SolverStyle(P, SamplingStyle(P)) == SolverStyle(paramP, SamplingStyle(P))
    @test sampling_grid(P,Pparam)≈sampling_grid(paramP, paramPparam)
    @test dimensions(dictionary(P,Pparam)) == dimensions(dictionary(paramP,paramPparam))
    @test dimensions(azdual_dict(P,Pparam))==dimensions(azdual_dict(paramP, paramPparam))
    @test AZ_A(P,Pparam)≈AZ_A(paramP, paramPparam)
    @test AZ_Zt(P,Pparam)≈AZ_Zt(paramP, paramPparam)
    @test AZ_A(P,Pparam)≈AZ_A(paramP, paramPparam)
    @test typeof(solver(P,Pparam))==typeof(solver(paramP,paramPparam))
end
