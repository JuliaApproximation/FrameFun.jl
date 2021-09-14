module test_suite_applications

using BasisFunctions, DomainSets, FrameFun, LinearAlgebra
using Test

@testset "FourierExtensionPlatform: DiscreteGramStyle" begin
    F = FourierExtensionPlatform(0.0..0.5)
    opts = (samplingstyle=FrameFun.DiscreteGramStyle(),)
    P = dictionary(F, 4)
    D = azdual_dict(F,4;opts...)
    M = measure(F,4;opts...)
    A = SynthesisOperator(P, M)
    Z = SynthesisOperator(D, M)
    Zt = Z'
    G1 = Zt*A
    G2 = BasisFunctions.default_mixedgram(D, P, M;warnslow=false)
    G3 = BasisFunctions.mixedgram(D, P, M)
    @test Matrix(G1) ≈ Matrix(G2)
    @test Matrix(G2) ≈ Matrix(G3)
end

@testset "FourierExtensionPlatform: GramStyle" begin
    F = FourierExtensionPlatform(0.0..0.5)
    opts = (samplingstyle=FrameFun.GramStyle(),)
    P = dictionary(F, 4)
    D = azdual_dict(F,4;opts...)
    M = measure(F,4;opts...)
    A = SynthesisOperator(P, M)
    Z = SynthesisOperator(D, M)
    @info "Two `Slow computation of Gram matrix entrywise.` warnings follow."
    Zt = Z'
    G1 = Zt*A
    G2 = BasisFunctions.default_mixedgram(D, P, M;warnslow=false)
    G3 = BasisFunctions.mixedgram(D, P, M)
    @test Matrix(G1) ≈ Matrix(G2)
    @test Matrix(G2) ≈ Matrix(G3)


    n = 101
    F = FourierExtensionPlatform(0.0..0.5)
    opts = (samplingstyle=FrameFun.GramStyle(),warnslow=false)
    A = AZ_A(F, n; opts...)
    Z = AZ_Z(F, n; opts...)
    MG = Z'*A
    @test sum(.0001 .< svdvals(MG) .< .9999) == 9
end

@testset "WeightedSumPlatform" begin
    # Simple platform construction test
    n = 12
    domain = 0.0..0.5
    P = FourierExtensionPlatform(domain)
    WP = WeightedSumPlatform(P,x->sqrt(x),x->1)
    f = x->sqrt(x)*(1-x)-exp(x)
    F = Fun(f, WP, (n,n), solverstyle=AZStyle())
    rgrid = randomgrid(domain, 200)
    abserror = sum(abs.(F.(rgrid)-f.(rgrid)))/length(rgrid)
    @test abserror<2e-7
end

end # module
