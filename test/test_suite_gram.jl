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
            b = FourierBasis{T}(n)
            d = interval(T(.25),T(.75))
            frame = extensionframe(b, d)
            g = Gram(Span(frame); reltol=tol,abstol=tol)
            m = matrix(g)
            @test norm(imag(m)) < tol
            m1 = real(m)

            I = [value(native_index(b, i)) for i in 1:n]
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
        B = ChebyshevBasis; n = 11; T = Float32;
        basis = instantiate(B, n, T)
        left(basis)
        for T in (Float32,Float64,), n in (10,11)
            tol = sqrt(eps(T))
            e = rand(T,n)
            for B in (ChebyshevBasis,LegendrePolynomials,FourierBasis,SineSeries,CosineSeries,BSplineTranslatesBasis,)
                basis = instantiate(B, n, T)
                domain = interval(left(basis),right(basis))
                frame = extensionframe(basis, domain)

                @test norm(Gram(Span(frame); abstol=tol, reltol=tol)*e - Gram(Span(basis); abstol=tol, reltol=tol)*e) <100*tol
                @test norm(DualGram(Span(frame); abstol=tol, reltol=tol)*e - DualGram(Span(basis); abstol=tol, reltol=tol)*e) <100*tol
                @test norm(MixedGram(Span(frame); abstol=tol, reltol=tol)*e - MixedGram(Span(basis); abstol=tol, reltol=tol)*e) <100*tol
            end
            for B in (ChebyshevBasis,FourierBasis,CosineSeries,BSplineTranslatesBasis,), oversampling in 1:4
                  basis = instantiate(B, n, T)
                  domain = interval(left(basis),right(basis))
                  frame = extensionframe(basis, domain)

                  @test norm(DiscreteGram(Span(frame); oversampling=oversampling)*e - DiscreteGram(Span(basis); oversampling=oversampling)*e) <100*tol
                  @test norm(DiscreteDualGram(Span(frame); oversampling=oversampling)*e - DiscreteDualGram(Span(basis); oversampling=oversampling)*e) <100*tol
                  @test norm(DiscreteMixedGram(Span(frame); oversampling=oversampling)*e - DiscreteMixedGram(Span(basis); oversampling=oversampling)*e) <100*tol
            end
        end
    end
end

function test_basis_oversampling()
    @testset "oversampling" begin
        for B in (ChebyshevBasis,CosineSeries,BSplineTranslatesBasis,), oversampling in 1:4, n in (10,11)
            basis = instantiate(B, n, Float64)
            domain = interval(left(basis),right(basis))
            @test BasisFunctions.basis_oversampling(extensionframe(basis, domain), oversampling)==oversampling
        end
        for B in (CosineSeries,BSplineTranslatesBasis,), oversampling in 1:4, n in (1000,1001)
            basis = instantiate(B, n, Float64)
            domain = interval(left(basis),right(basis))/2
            # println(FrameFun.basis_oversampling(domain, basis, oversampling), " ", 2oversampling)
            @test abs(BasisFunctions.basis_oversampling(extensionframe(basis, domain), oversampling)-2oversampling) < .01
        end
    end
end

function test_discrete_gram()

    @testset "Testing discrete dual gram en mixed gram with oversampling" begin
        for B in [ChebyshevBasis, FourierBasis, BSplineTranslatesBasis]
            b = instantiate(B, 1, Float64)
            d = interval(left(b), right(b))/2
            for T in [Float64, BigFloat]
                for n in [10,11]
                    e = map(T, rand(n))
                    b = instantiate(B,n,T)
                    for os in 1:4

                        frame = extensionframe(b,d)
                        Gomega = DiscreteGram(Span(frame); oversampling=os)
                        Eomega = evaluation_operator(Span(frame); oversampling=os)
                        N = BasisFunctions.discrete_gram_scaling(frame, os)

                        basis_os = BasisFunctions.basis_oversampling(frame,os)

                        @test (Eomega'Eomega)*e/N≈Gomega*e

                        GT = DiscreteDualGram(Span(b); oversampling=basis_os)

                        ETomega = Eomega*GT

                        GTomega = GT*Gomega*GT

                        @test (ETomega'ETomega)*e/N≈GTomega*e

                        GMomega = GT*Gomega

                        @test (ETomega'Eomega)*e/N≈GMomega*e

                        @test DiscreteGram(Span(frame);oversampling=os)*e≈Gomega*e
                        @test DiscreteDualGram(Span(frame); oversampling=os)*e≈GTomega*e
                        @test DiscreteMixedGram(Span(frame); oversampling=os)*e≈GMomega*e
                    end
                end
            end
        end
    end
end

function test_connection_restriction_extension_discretegram()
  @testset "Testing connection extension with discrete gram" begin
    for T in (Float64, BigFloat), (b,os) in [(BSplineTranslatesBasis(10, 1, T),1), (BSplineTranslatesBasis(10, 2, T),1), (FourierBasis{T}(10),1.1), (ChebyshevBasis{T}(20),.66)]
      d = interval(left(b),right(b))/2
      frame = extensionframe(b, d)
      bspan = Span(b)
      framespan = Span(frame)
      N = BasisFunctions.discrete_gram_scaling(b, os)

       # Uses extension times two next, so works only for os=1, for bsplines, ≈ 1.1 for fourier series,...
      @assert 2≈BasisFunctions.basis_oversampling(frame, os)
      G = DiscreteGram(framespan; oversampling=os)

      # check whether the previous calculation is the same as extension by 2.
      b_large = extend(b)
      b_largespan = Span(b_large)
      time_basisspan = gridspace(b_largespan)
      r_time_basisspan = gridspace(b_largespan,FrameFun.subgrid(grid(b_large), d))


      E = extension_operator(bspan, b_largespan)
      A = evaluation_operator(b_largespan)
      R = restriction_operator(time_basisspan, r_time_basisspan)
      A_Omega = R*A*E
      matrix(A_Omega)
      e = map(T, rand(size(A_Omega,2)))
      @assert A_Omega*e ≈ evaluation_operator(framespan, oversampling=os)*e
      G_test = (1/T(N))*A_Omega'A_Omega

      e = map(T, rand(size(G,2)))
      @test 1+maximum(abs.((G- G_test)*e))≈1

      GD = DiscreteDualGram(framespan; oversampling=os)

      Ad = discrete_dual_evaluation_operator(b_largespan)
      Ad_Omega = R*discrete_dual_evaluation_operator(bspan; oversampling=oversampling=BasisFunctions.basis_oversampling(frame,os))
      GD_test = (1/T(N))*Ad_Omega'Ad_Omega

      @test 1+maximum(abs.((GD- GD_test)*e))≈1.

      GM = DiscreteMixedGram(framespan; oversampling=os)
      GM_test = (1/T(N))*Ad_Omega'A_Omega

      @test 1+maximum(abs.((GM- GM_test)*e))≈1
    end
  end
  for T in (Float64,BigFloat), B in (BSplineTranslatesBasis, FourierBasis, ChebyshevBasis,), n in (10,11), os in 1:4
    b = instantiate(B, n, T)
    bspan = Span(b)
    d = interval(left(b),right(b))/2
    frame = extensionframe(b, d)
    framespan = Span(frame)
    N = BasisFunctions.discrete_gram_scaling(frame, os)

    G = DiscreteGram(framespan; oversampling=os)

    r_grid = BasisFunctions.oversampled_grid(frame, os)
    t_grid = BasisFunctions.oversampled_grid(b, BasisFunctions.basis_oversampling(frame, os))
    time_basisspan = gridspace(bspan, t_grid)
    r_time_basisspan = gridspace(framespan, r_grid)

    A = evaluation_operator(bspan; oversampling=BasisFunctions.basis_oversampling(frame, os))
    R = restriction_operator(time_basisspan, r_time_basisspan)
    Af = evaluation_operator(framespan; oversampling=os)
    e = map(T, rand(size(A,2)))
    @test (R*A)*e ≈ Af*e
    A_Omega = Af
    G_test = (1/T(N))*A_Omega'A_Omega

    e = map(T, rand(size(G,2)))
    @test 1+maximum(abs.((G- G_test)*e))≈1

    GD = DiscreteDualGram(framespan; oversampling=os)

    Ad = discrete_dual_evaluation_operator(bspan; oversampling=oversampling=BasisFunctions.basis_oversampling(frame,os))
    Adf = discrete_dual_evaluation_operator(framespan, oversampling=os)
    e = map(T, rand(size(Ad,2)))
    @test (R*Ad)*e ≈ Adf*e

    Ad_Omega = Adf
    GD_test = (1/T(N))*Ad_Omega'Ad_Omega

    @test 1+maximum(abs.((GD- GD_test)*e))≈1.

    GM = DiscreteMixedGram(framespan; oversampling=os)
    GM_test = (1/T(N))*Ad_Omega'A_Omega

    @test 1+maximum(abs.((GM- GM_test)*e))≈1
  end
end

function test_general_gram()
    T = Float64
    tol = max(sqrt(eps(T)), 1e-10)
    for method in (Gram, DualGram, MixedGram), B in (FourierBasis{T}(11), BSplineTranslatesBasis(5, 1,T))
        D = interval(left(B), right(B))
        frame = extensionframe(B,D)
        GBB = method(Span(frame),Span(frame); abstol=tol, reltol=tol)
        GB = method(Span(B); abstol=tol, reltol=tol)

        e = rand(length(frame))
        @test norm(GBB*e - GB*e) <= 1000*tol
    end
end

delimit("Gram")
basis_test()
prolate_test()
test_basis_oversampling()
test_discrete_gram()
test_connection_restriction_extension_discretegram()
test_general_gram()

end
