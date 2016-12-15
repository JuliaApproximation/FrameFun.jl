module test_suite


using BasisFunctions
using FrameFun
using Base.Test
FE = FrameFun
BA = BasisFunctions

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

function test_function_space()
  bboxes = (Interval(), Interval(),
      Interval(0,1), Interval(), Interval(),
      Interval(0,1)⊗Interval(),Interval(-2,1),
      Interval(0,1)⊗Interval(0,1)⊗Interval(0,1), Interval(0,1))
  bases = (FourierBasis(121,-1,1), FourierBasis(121,-1,1),
      FourierBasis(121), FourierBasis(121,-1,1), ChebyshevBasis(121),
      FourierBasis(121)⊗ChebyshevBasis(121), FourierBasis(121,-2.,1.)⊕rescale(ChebyshevBasis(121),-2.,1.),
      BA.tensorproduct(FourierBasis(121),3), BA.multiset(FourierBasis(121),FourierBasis(121),FourierBasis(121)))
  @testset "Space = $(name(space)) " for (i,space) in enumerate([
      FE.FunctionSpace(FourierBasis(121,-1,1)),FE. FunctionSpace(FourierBasis(121,-1,1), FE.BBox(-1,1)),
      FourierSpace(), FourierSpace(-1,1), ChebyshevSpace(),
      FourierSpace()⊗ChebyshevSpace(),FourierSpace(-2,1)⊕ChebyshevSpace(-2,1),FE.tensorproduct(FourierSpace(),3), FE.add(FourierSpace(),3)])
      @test left(bboxes[i])==left(boundingbox(space))
      @test right(bboxes[i])==right(boundingbox(space))
      @test FunctionSet(space, 121) == bases[i]
  end
  @testset "Util functions" begin
    for n in 1:4
      S = FE.tensorproduct(FourierSpace(),n)
      @test ndims(typeof(S)) == n
      @test ndims(S) == n
      @test numtype(S) == Float64
      @test eltype(S) == Complex128
    end
    @test eltype(promote_eltype(ChebyshevSpace(),Complex128)) == Complex128
    @test eltype(promote_eltype(ChebyshevSpace(),Float32)) == Float64
    @test eltype(promote_eltype(ChebyshevSpace(),Float64)) == Float64
    @test promote(ChebyshevSpace(), ChebyshevSpace()) == (ChebyshevSpace(), ChebyshevSpace())
    @test promote(ChebyshevSpace(), FourierSpace()) == (promote_eltype(ChebyshevSpace(),Complex128), FourierSpace())
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
        resnew = FE.residual(f,F)
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
    for tol in 10.0.^(-4.:-4.:-16.)
        F = mth(f, S, D, tol=tol, max_logn_coefs=max_logn_coefs)
        @test ( (FE.maxerror(f,F) < tol*100) || (length(set(F))>= 2^max_logn_coefs))
    end
    F = mth(x->cos(4π*x), S, D)
    @test FE.residual(x->cos(4π*x),F) < 1e-13
  end
end

function test_funs2d()
  tol = 1e-8
  D = Cube([0,0],[0.5,0.5])
  FF = FourierSpace()⊗FourierSpace()
  f(x,y) = exp(y*2*x)
  FC = FunConstructor(FF, D)
  F0 = FC(f,tol=tol, max_logn_coefs=12)
  @testset "test funs in 2D" begin @test FE.maxerror(f,F0) < 100*tol end
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

delimit("Adaptivity")
test_function_space()
test_residual()
test_funs()
test_funs2d()
test_extra_functionality()
end
