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
      for B in (ChebyshevBasis,FourierBasis,CosineSeries,BSplineTranslatesBasis,), oversampling in 1:4
        basis = instantiate(B, n, T)
        domain = Interval(left(basis),right(basis))
        frame = extensionframe(basis, domain)

        @test norm(DiscreteGram(frame; oversampling=oversampling)*e - DiscreteGram(basis; oversampling=oversampling)*e) <100*tol
        @test norm(DiscreteDualGram(frame; oversampling=oversampling)*e - DiscreteDualGram(basis; oversampling=oversampling)*e) <100*tol
        @test norm(DiscreteMixedGram(frame; oversampling=oversampling)*e - DiscreteMixedGram(basis; oversampling=oversampling)*e) <100*tol
      end
    end
  end
end

function test_basis_oversampling()
  @testset "oversampling" begin
    for B in (ChebyshevBasis,CosineSeries,BSplineTranslatesBasis,), oversampling in 1:4, n in (10,11)
      basis = instantiate(B, n, Float64)
      domain = Interval(left(basis),right(basis))
      @test BasisFunctions.basis_oversampling(extensionframe(basis, domain), oversampling)==oversampling
    end
    for B in (CosineSeries,BSplineTranslatesBasis,), oversampling in 1:4, n in (1000,1001)
      basis = instantiate(B, n, Float64)
      domain = Interval(left(basis),right(basis))/2
      # println(FrameFun.basis_oversampling(domain, basis, oversampling), " ", 2oversampling)
      @test abs(BasisFunctions.basis_oversampling(extensionframe(basis, domain), oversampling)-2oversampling) < .01
    end
  end
end

function test_discrete_gram()
  @testset "Testing discrete dual gram en mixed gram with oversampling" begin
    for T in [Float64, BigFloat]
      for n in [10,11], os in 1:4, B in [ChebyshevBasis, FourierBasis, BSplineTranslatesBasis]
        e = map(T, rand(n))
        b = instantiate(B,n,T)
        d = Interval(left(b), right(b))/2
        frame = extensionframe(b,d)
        Gomega = DiscreteGram(frame; oversampling=os)
        Eomega = evaluation_operator(frame; oversampling=os)
        N = length(BasisFunctions.oversampled_grid(frame, os))

        basis_os = BasisFunctions.basis_oversampling(frame,os)


        @test (Eomega'Eomega)*e/N≈Gomega*e

        GT = DiscreteDualGram(b; oversampling=basis_os)

        ETomega = Eomega*GT

        GTomega = GT*Gomega*GT

        @test (ETomega'ETomega)*e/N≈GTomega*e

        GMomega = GT*Gomega

        @test (ETomega'Eomega)*e/N≈GMomega*e

        @test DiscreteGram(frame;oversampling=os)*e≈Gomega*e
        @test DiscreteDualGram(frame; oversampling=os)*e≈GTomega*e
        @test DiscreteMixedGram(frame; oversampling=os)*e≈GMomega*e
      end
    end
  end
end

function test_connection_restriction_extension_discretegram()
  @testset "Testing connection extension with discrete gram" begin
    for T in (Float64, BigFloat), (b,os) in [(BSplineTranslatesBasis(10, 1, T),1), (BSplineTranslatesBasis(10, 2, T),1), (FourierBasis(10,T),1.1), (ChebyshevBasis(20,T),.66)]
      d = Interval(left(b),right(b))/2
      frame = extensionframe(b, d)

       # Uses extension times two next, so works only for os=1, for bsplines, ≈ 1.1 for fourier series,...
      @assert 2≈BasisFunctions.basis_oversampling(frame, os)
      G = DiscreteGram(frame; oversampling=os)

      # check whether the previous calculation is the same as extension by 2.
      b_large = extend(b)
      time_basis = gridspace(b_large)
      r_time_basis = gridspace(b_large,FrameFun.subgrid(grid(b_large), d))



      E = extension_operator(b, b_large)
      A = evaluation_operator(b_large)
      R = restriction_operator(time_basis, r_time_basis)
      A_Omega = R*A*E
      matrix(A_Omega)
      e = map(T, rand(size(A_Omega,2)))
      @assert A_Omega*e ≈ evaluation_operator(frame, oversampling=os)*e
      G_test = (1/T(length(r_time_basis)))*A_Omega'A_Omega

      e = map(T, rand(size(G,2)))
      @test 1+maximum(abs((G- G_test)*e))≈1

      GD = DiscreteDualGram(frame; oversampling=os)

      Ad = discrete_dual_evaluation_operator(b_large)
      # Ad_Omega = R*Ad*E
      Ad_Omega = R*discrete_dual_evaluation_operator(b; oversampling=oversampling=BasisFunctions.basis_oversampling(frame,os))
      GD_test = (1/T(length(r_time_basis)))*Ad_Omega'Ad_Omega

      @test 1+maximum(abs((GD- GD_test)*e))≈1.

      GM = DiscreteMixedGram(frame; oversampling=os)
      GM_test = (1/T(length(r_time_basis)))*Ad_Omega'A_Omega

      @test 1+maximum(abs((GM- GM_test)*e))≈1
    end
  end
  for T in (Float64,BigFloat), B in (BSplineTranslatesBasis, FourierBasis, ChebyshevBasis,), n in (10,11), os in 1:4
    b = instantiate(B, n, T)
    d = Interval(left(b),right(b))/2
    frame = extensionframe(b, d)

    G = DiscreteGram(frame; oversampling=os)

    r_grid = BasisFunctions.oversampled_grid(frame, os)
    t_grid = BasisFunctions.oversampled_grid(b, BasisFunctions.basis_oversampling(frame, os))
    time_basis = gridspace(b, t_grid)
    r_time_basis = gridspace(frame, r_grid)

    A = evaluation_operator(b; oversampling=BasisFunctions.basis_oversampling(frame, os))
    R = restriction_operator(time_basis, r_time_basis)
    Af = evaluation_operator(frame; oversampling=os)
    e = map(T, rand(size(A,2)))
    @test (R*A)*e ≈ Af*e
    A_Omega = Af

    G_test = (1/T(length(r_time_basis)))*A_Omega'A_Omega

    e = map(T, rand(size(G,2)))
    @test 1+maximum(abs((G- G_test)*e))≈1

    GD = DiscreteDualGram(frame; oversampling=os)

    Ad = discrete_dual_evaluation_operator(b; oversampling=oversampling=BasisFunctions.basis_oversampling(frame,os))
    Adf = discrete_dual_evaluation_operator(frame, oversampling=os)
    e = map(T, rand(size(Ad,2)))
    @test (R*Ad)*e ≈ Adf*e

    Ad_Omega = Adf
    GD_test = (1/T(length(r_time_basis)))*Ad_Omega'Ad_Omega

    @test 1+maximum(abs((GD- GD_test)*e))≈1.

    GM = DiscreteMixedGram(frame; oversampling=os)
    GM_test = (1/T(length(r_time_basis)))*Ad_Omega'A_Omega

    @test 1+maximum(abs((GM- GM_test)*e))≈1
  end
end

delimit("Gram")
basis_test()
prolate_test()
test_basis_oversampling()
test_discrete_gram()
test_connection_restriction_extension_discretegram()

end
