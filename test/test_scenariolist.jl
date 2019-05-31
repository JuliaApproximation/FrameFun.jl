using BasisFunctions, LinearAlgebra, DomainSets, Grids, Test, StaticArrays, FrameFun
@testset begin
    B = Fourier(11,-1,1)⊗Fourier(11,-1,1)
    Dom = disk(0.8)
    @test support(dictionary(∂x(random_expansion(extensionframe(B, Dom)))))≈Dom
end
