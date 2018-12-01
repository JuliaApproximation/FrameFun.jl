module test_suite

using DomainSets, BasisFunctions, FrameFun

using Test, LinearAlgebra, Random
Random.seed!(1234)

FE = FrameFun
BA = BasisFunctions

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

function test_function_space()
    bboxes = (Interval(-1.0, 1.0),
        Interval(-1.0, 1.0),
        Interval(0.0, 1.0),
        Interval(-1.0, 1.0),
        Interval(-1.0, 1.0),
        Interval(0.0, 1.0)×Interval(-1.0, 1.0),
        Interval(-2.0, 1.0),
        Interval(0.0, 1.0)×Interval(0.0, 1.0)×Interval(0.0, 1.0),
        Interval(0.0, 1.0))
    bases = (FourierBasis(64,-1.0, 1.0),
        FourierBasis(64,-1.0, 1.0),
        FourierBasis(64),
        FourierBasis(64,-1.0, 1.0),
        ChebyshevBasis(64),
        FourierBasis(8)⊗ChebyshevBasis(8),
        FourierBasis(32,-2.,1.)⊕rescale(ChebyshevBasis(32),-2.,1.),
        BA.tensorproduct(FourierBasis(4),3),
        BA.multidict(FourierBasis(32),FourierBasis(32)))
    @testset "Space = $(name(space)) " for (i,space) in enumerate([
        FE.FunctionSpace(FourierBasis(64,-1.0, 1.0)),
        FE.FunctionSpace(FourierBasis(64,-1.0, 1.0), Interval(-1.0, 1.0)),
        FourierSpace(),
        FourierSpace(-1.0, 1.0),
        ChebyshevSpace(),
        FourierSpace() ⊗ ChebyshevSpace(),
        FourierSpace(-2,1) ⊕ ChebyshevSpace(-2,1),
        FE.tensorproduct(FourierSpace(),3),
        FE.add(FourierSpace(),2)])
        # @test left(bboxes[i])==left(boundingbox(space))
        # @test right(bboxes[i])==right(boundingbox(space))
        @test Dictionary(space, 64) == bases[i]
    end
    @testset "Util functions" begin
    for n in 1:4
        S = FE.tensorproduct(FourierSpace(),n)
        @test dimension(S) == n
        @test codomaintype(S) == ComplexF64
    end
    @test domaintype(promote_domaintype(ChebyshevSpace(),ComplexF64)) == domaintype(promote_domaintype(ChebyshevBasis(11),ComplexF64))
    @test domaintype(promote_domaintype(ChebyshevSpace(),Float32)) == domaintype(promote_domaintype(ChebyshevBasis(11),Float32))
    @test domaintype(promote_domaintype(ChebyshevSpace(),Float64)) == domaintype(promote_domaintype(ChebyshevBasis(11),Float64))
    end
end

function test_residual()
    @testset "Residual for basis $(basis)" for basis in (FourierBasis, ChebyshevBasis)
        domain = Interval(-1.0, 1.0)/2
        f = x->cos(20x)
        res = Inf
        for n in 2 .^ (3:5)
            S = rescale(instantiate(basis,n), -1.0, 1.0)
            F = Fun(S, domain, f)
            resnew = FE.residual(f,F)
            @test resnew < res
            res = resnew
        end
    end
end

function test_funs()
    tests = ("fun_simple", "fun_optimal_N", "fun_greedy")
    S = FourierBasis(0, -1.0, 1.0)
    D = Interval(-1.0, 1.0)/2
    f = x->x
    max_logn_coefs = 8
    @testset "$(tests[i]) tests" for (i,mth) in enumerate([FE.fun_simple, FE.fun_optimal_N, FE.fun_greedy])
        for tol in 10.0.^(-4.:-4.:-12.)
            F = mth(f, S, D; tol=tol, max_logn_coefs=max_logn_coefs, adaptive_verbose=false)
            @test (FE.maxerror(f,F) < tol*100) || (length(dictionary(F))>= 2^max_logn_coefs)
        end
        if mth != FE.fun_greedy
            for tol in 10.0.^(-4.:-4.:-16.)
                F = mth(f, S, D; tol=tol, max_logn_coefs=max_logn_coefs, abscoef=3, adaptive_verbose=false, random_test=true)
                @test (FE.maxerror(f,F) < tol*100) || (length(dictionary(F))>= 2^max_logn_coefs)
            end
        end
        F = mth(x->cos(4π*x), S, D)
        @test residual(x->cos(4π*x), F, residualtype = :l2) < 1e-11
    end
end

function test_funs2d()
    tol = 1e-8
    D = cube([0.0, 0.0], [0.5, 0.5])
    FF = FourierSpace()⊗FourierSpace()
    f(x,y) = exp(y*2*x)
    FC = FunConstructor(FF, D)
    F0 = FC(f,tol=tol, max_logn_coefs=12)
    @testset "test funs in 2D" begin @test FE.maxerror(f,F0) < 100*tol end
end

function test_extra_functionality()
    F = FourierSpace(-1.0, 1.0)
    D = Interval(-0.5, 0.5)
    FC = FunConstructor(F, D)
    x = FC(identity)
    tol = 1e-9
    @testset "Extra Functionality" begin
        f4(x) = sin(cos(x))
        F4 = sin(cos(x,tol=tol))
        @test(FE.maxerror(f4,F4) < 10*tol)
        f5(x) = exp(cos(20x))
        F5 = exp(cos(20x, afun=FE.fun_simple, tol=tol))
        @test(FE.maxerror(f5,F5) < 10*tol)
    end
end

delimit("Adaptivity")

test_function_space()
test_residual()
test_funs()
test_funs2d()
test_extra_functionality()

end
