module test_suite_adaptive

using BasisFunctions
using FrameFuns
using Base.Test
FE = FrameFuns
BA = BasisFunctions

function test_function_space()
  set = FourierBasis(121,-1,1)
  bboxes = (Interval(), Interval(),
      Interval(0,1), Interval(), Interval(),
      Interval(0,1)⊗Interval(),
      Interval(-2,1))
  bases = (set, set,
      FourierBasis(121), FourierBasis(121,-1,1), ChebyshevBasis(121),
      FourierBasis(121)⊗ChebyshevBasis(121), FourierBasis(121,-2.,1.)⊕rescale(ChebyshevBasis(121),-2.,1.))
  @testset "Space = $(name(space)) " for (i,space) in enumerate([
      FE.FunctionSpace(set),FE. FunctionSpace(set, FE.BBox(-1,1)),
      FourierSpace(), FourierSpace(-1,1), ChebyshevSpace(),
      FourierSpace()⊗ChebyshevSpace(), FourierSpace(-2,0)⊕ChebyshevSpace()])
    @test left(bboxes[i])==left(boundingbox(space))
    @test right(bboxes[i])==right(boundingbox(space))
    @test FunctionSet(space, 121) == bases[i]
  end
end

function test_residual()
  @testset "Residual for basis $(basis)" for basis in (FourierBasis, ChebyshevBasis)
    D = Interval()/2
    f = x->cos(20x)
    res = Inf
    for n in 2.^(3:6)
        S = rescale(instantiate(basis,n), -1,1)
        F = Fun(f, S, D)
        resnew = FE.residual(F,f)
        @test  resnew < res
        res = resnew
    end
  end
end

function test_funs()
  tests = ("fun_simple", "fun_optimal_N", "fun_greedy")
  S = FourierBasis(0,-1,1)
  D = Interval()/2
  f = x->x
  max_logn_coefs = 8
  @testset "$(tests[i]) tests" for (i,mth) in enumerate([FE.fun_simple, FE.fun_optimal_N, FE.fun_greedy])
    for tol in 10.0.^(-6.:-2.:-14.)
        F = mth(f, S, D, tol=tol, max_logn_coefs=max_logn_coefs)
        @test ( (FE.maxerror(f,F) < tol*100) || (length(set(F))>= 2^max_logn_coefs))
    end
    F = mth(x->cos(4π*x), S, D)
    @test FE.residual(x->cos(4π*x),F) < 1e-13
  end
end

function test_extra_functionality()
  F = FourierSpace(-1,1)
  D = Interval(-0.5,0.5)
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
println()
println("############")
println("# Adaptivity")
println("############")
test_function_space()
test_residual()
test_funs()
test_extra_functionality()
end
