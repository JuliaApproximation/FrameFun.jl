module test_suite

using BasisFunctions
using FrameFun
using Base.Test
FE = FrameFun
BA = BasisFunctions

function prolate_test()
  @testset "Fourier frame on half, symmetric domain (prolates)" begin
    for T in (Float32, Float64,), n in (11,13,)
      tol = sqrt(eps(T))
      b = FourierBasis(n,T)
      d = Interval(T(.25),T(.75))
      frame = extensionframe(b, d)
      g = Gram(frame; reltol=tol,abstol=tol)
      m = matrix(g)
      @test norm(imag(m)) < tol
      m1 = real(m)

      I = [native_index(b, i).index for i in 1:n]
      m2 = real([k==l ? .5: exp(1im*pi*(k-l))/(pi*(k-l))*sin(2pi*(1/2-T(.25))*(k-l)) for k in I, l in I])
      @test norm(m1-m2) < tol
    end
  end
end

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

function basis_test()
  @testset "Frame on entire domain (i.e, basis problem)" begin
    for T in (Float32,Float64,), n in (10,11)
      tol = sqrt(eps(T))
      e = rand(T,n)
      for B in (ChebyshevBasis,LegendreBasis,FourierBasis,SineSeries,CosineSeries,BSplineTranslatesBasis,)
        basis = instantiate(B, n, T)
        domain = Interval(left(basis),right(basis))
        frame = extensionframe(basis, domain)

        @test norm(Gram(frame; abstol=tol, reltol=tol)*e - Gram(basis; abstol=tol, reltol=tol)*e) <100*tol
        @test norm(DualGram(frame; abstol=tol, reltol=tol)*e - DualGram(basis; abstol=tol, reltol=tol)*e) <100*tol
        @test norm(MixedGram(frame; abstol=tol, reltol=tol)*e - MixedGram(basis; abstol=tol, reltol=tol)*e) <100*tol
      end
      for B in (ChebyshevBasis,FourierBasis,CosineSeries,BSplineTranslatesBasis,)
        basis = instantiate(B, n, T)
        domain = Interval(left(basis),right(basis))
        frame = extensionframe(basis, domain)

        @test norm(DiscreteGram(frame)*e - DiscreteGram(basis)*e) <100*tol
        @test norm(DiscreteDualGram(frame)*e - DiscreteDualGram(basis)*e) <100*tol
        @test norm(DiscreteMixedGram(frame)*e - DiscreteMixedGram(basis)*e) <100*tol
      end
    end
  end
end

delimit("Gram")
basis_test()
prolate_test()


end
