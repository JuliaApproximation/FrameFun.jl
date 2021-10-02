using Test, LinearAlgebra, BasisFunctions, FrameFun

dict1 = Fourier(10)
dict2 = Fourier(10)^2
plat1 = platform(dict1)
plat2 = platform(dict2)
par1 = 10
par2 = (10,10)

all_dicts = (dict1, dict2)

@testset "samplingparameter" begin
    # adaptive approximations has no samplingparameter
    for (a,sol) in zip(all_dicts ,(10,(10,10),10,(10,10)))
        @test samplingparameter(a)==sol
    end
    @test samplingparameter(plat1, par1) == 10
    @test samplingparameter(plat2, par2) == (10,10)
end

@testset "interpolation_grid" begin
    @test interpolation_grid(plat1, par1) == interpolation_grid(dict1)
    @test interpolation_grid(plat2, par2) == interpolation_grid(dict2)

    ap = approximationproblem(dict1)
    @test !FrameFun.ap_hasproperty(ap, interpolation_grid)
    @test interpolation_grid(ap) == interpolation_grid(dict1)
    @test FrameFun.ap_hasproperty(ap, interpolation_grid)
    @test FrameFun.cache(ap, interpolation_grid) == interpolation_grid(dict1)
end

@testset "oversampling_grid" begin
    @test oversampling_grid(dict1) == FourierGrid(10)
    @test oversampling_grid(dict1, 11) == FourierGrid(11)
    @test oversampling_grid(dict2) == FourierGrid(10)^2
    @test oversampling_grid(dict2, (11,11)) == FourierGrid(11)^2
    @test oversampling_grid(plat1, par1) == FourierGrid(10)

    ap = FrameFun.approximationproblem(dict1, 15)
    @test !FrameFun.ap_hasproperty(ap, oversampling_grid)
    @test oversampling_grid(ap) == FourierGrid(15)
    @test FrameFun.ap_hasproperty(ap, oversampling_grid)
    @test FrameFun.cache(ap, oversampling_grid) == FourierGrid(15)
end

@testset "sampling_grid" begin
    @test sampling_grid(dict1) == FourierGrid(10)
    @test sampling_grid(dict2) == FourierGrid(10)^2
    @test sampling_grid(plat1, par1) == FourierGrid(10)
    @test sampling_grid(plat2, par2) == FourierGrid(10)^2
end

@testset "platform_grid" begin
    @test platform_grid(dict1; samplingstyle=GridStyle(), grid=FourierGrid(1)) == FourierGrid(1)
    @test platform_grid(dict2; samplingstyle=GridStyle(), grid=FourierGrid(1)^2) == FourierGrid(1)^2
    @test platform_grid(plat1, par1, grid=FourierGrid(1)) == FourierGrid(1)
    @test sampling_grid(dict2; samplingstyle=GridStyle(),grid=FourierGrid(1)^2)==FourierGrid(1)^2
end

@testset "sampling_operator" begin
    op1 = sampling_operator(dict1)
    op2 = sampling_operator(dict2)
    op3 = sampling_operator(plat1, par1)
    op4 = sampling_operator(plat2, par2)

    @test grid(dest(op1)) == FourierGrid(10)
    @test grid(dest(op3)) == FourierGrid(10)

    f2 = (x,y)->exp(x*y)
    @test tensorproduct(op1,op1)*f2≈op2*f2
    @test tensorproduct(op3,op3)*f2≈op4*f2
end

@testset "discretemeasure" begin
    @test discretemeasure(dict1) == discretemeasure(sampling_grid(dict1))
    @test discretemeasure(dict2)== discretemeasure(sampling_grid(dict2))
    @test discretemeasure(plat1, par1)== discretemeasure(sampling_grid(plat1, par1))
    @test discretemeasure(plat2, par2)== discretemeasure(sampling_grid(plat2, par2))
end

@testset "measure" begin
    @test measure(platform(dict1)) == measure(dict1)
    @test measure(platform(dict2)) == measure(dict2)
end

@testset "azdual" begin
    @test operator(azdual(dict1))≈operator(azdual(plat1, par1))
    op = azdual(dict2)
    @test op isa TensorProductDict
    @test ≈(operator.(components(op))...)
    op = azdual(plat2, par2)
    @test op isa TensorProductDict
    @test ≈(operator.(components(op))...)
    @test operator(component(azdual(dict2),1))≈operator(azdual(dict1))

    @test operator(azdual(dict1;samplingstyle=GramStyle())) ≈
        operator(component(azdual(dict2;samplingstyle=ProductSamplingStyle(GramStyle(),GramStyle())),1))
    @test operator(azdual(dict1;samplingstyle=GramStyle())) ≈
        operator(component(azdual(plat2, par2;samplingstyle=ProductSamplingStyle(GramStyle(),GramStyle())),1))
    @test operator(azdual(dict1;samplingstyle=GramStyle())) ≈
        operator(component(azdual(dict2;samplingstyle=ProductSamplingStyle(GramStyle(),GramStyle())),1))
end

@testset "discretization" begin
    op1 = discretization(plat1, par1)
    op2 = discretization(plat2, par2)
    op3 = discretization(dict1)
    op4 = discretization(dict2)

    @test op2 isa TensorProductOperator
    @test op4 isa TensorProductOperator

    @test op1 ≈ op3
    @test op2 ≈ op4

    f1 = exp
    f2 = (x,y)->exp(x*y)
    op3, b3 = full_discretization(f1, plat1, par1)
    op4, b4 = full_discretization(f2, plat2, par2)
    op5, b5 = full_discretization(f1, dict1)
    op6, b6 = full_discretization(f2, dict2)

    @test op4 isa TensorProductOperator
    @test op6 isa TensorProductOperator

    @test op3 ≈ op5
    @test op4 ≈ op6

    @test b3≈b5
    @test b4≈b6
end

@testset "dualdiscretization" begin
    op1 = dualdiscretization(plat1, par1)
    op2 = dualdiscretization(plat2, par2)
    op3 = dualdiscretization(dict1)
    op4 = dualdiscretization(dict2)

    @test op2 isa TensorProductOperator
    @test op4 isa TensorProductOperator

    @test op1 ≈ op3
    @test op2 ≈ op4
end

@testset "solver" begin
    op1 = solver(dict1)
    op2 = solver(dict2)
    op3 = solver(plat1, par1)
    op4 = solver(plat2, par2)

    @test op2 isa TensorProductOperator
    @test op4 isa TensorProductOperator

    @test op1 ≈ op3
    @test op2 ≈ op4

    @test component(op2,1) ≈ op1 ≈ component(op2,2)

    for STYLE in (InverseStyle, DirectStyle, DualStyle, IterativeStyle,AZStyle,AZSmoothStyle)
        op1 = solver(dict1;solverstyle=STYLE())
        op2 = solver(dict2;solverstyle=ProductSolverStyle(STYLE(),STYLE()))
        op3 = solver(plat1, par1;solverstyle=STYLE())
        op4 = solver(plat2, par2;solverstyle=ProductSolverStyle(STYLE(),STYLE()))

        @test op2 isa TensorProductOperator
        @test op4 isa TensorProductOperator

        @test op1 ≈ op3
        @test op2 ≈ op4

        @test component(op2,1) ≈ op1 ≈ component(op2,2)
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
    P = ExtensionFramePlatform(WeightedSumPlatform(platform(ChebyshevT(10)^2), (x,y)->1.,
                (x,y)->sqrt(x^2+y^2)), .9 .* UnitDisk())
    path = HierarchyPath(p,ProductPath(p,p))
    paramP = parametrizedplatform(P, path)

    paramPparam=45
    Pparam = path[paramPparam]

    @test SamplingStyle(P)==SamplingStyle(paramP)
    @test SolverStyle(P, SamplingStyle(P)) == SolverStyle(paramP, SamplingStyle(P))
    @test sampling_grid(P,Pparam)≈sampling_grid(paramP, paramPparam)
    @test dimensions(dictionary(P,Pparam)) == dimensions(dictionary(paramP,paramPparam))
    @test dimensions(azdual(P,Pparam))==dimensions(azdual(paramP, paramPparam))
    @test AZ_A(P,Pparam)≈AZ_A(paramP, paramPparam)
    @test AZ_Zt(P,Pparam)≈AZ_Zt(paramP, paramPparam)
    @test AZ_A(P,Pparam)≈AZ_A(paramP, paramPparam)
    @test typeof(solver(P,Pparam))==typeof(solver(paramP,paramPparam))
end
