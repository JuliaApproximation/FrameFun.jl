using Test, FrameFun, DomainSets, BasisFunctions

@testset "FrameFunPlatforms" begin
    plat1 = FourierPlatform()
    plat2 = FourierExtensionPlatform(0.0..0.5)
    plat3 = ChebyshevPlatform()
    plat4 = OPSExtensionFramePlatform(ChebyshevU(10),0.0..0.5)

    measure(plat1) == FourierWeight()
    measure(plat2) == submeasure(FourierWeight(),0.0..0.5)
    measure(plat3) == ChebyshevWeight()
    measure(plat4) == submeasure(ChebyshevWeight(),0.0..0.5)

    @test SamplingStyle(plat1) == InterpolationStyle()
    @test SamplingStyle(plat2) == OversamplingStyle()
    @test SamplingStyle(plat3) == InterpolationStyle()
    @test SamplingStyle(plat4) == OversamplingStyle()

    @test DictionaryStyle(plat1) == BasisStyle()
    @test DictionaryStyle(plat2) == FrameStyle()
    @test DictionaryStyle(plat3) == BasisStyle()
    @test DictionaryStyle(plat4) == FrameStyle()

    @test SolverStyle(plat1,InterpolationStyle()) == InverseStyle()
    @test SolverStyle(plat2,OversamplingStyle()) == AZStyle()
    @test SolverStyle(plat3,InterpolationStyle()) == InverseStyle()
    @test SolverStyle(plat4,OversamplingStyle()) == AZStyle()

    @test dictionary(plat1,10) == Fourier(10)
    @test iscompatible(dictionary(plat2,10), extensionframe(Fourier(10),0.0..0.5))
    @test dictionary(plat3,10) == ChebyshevT(10)
    @test iscompatible(dictionary(plat4,10), extensionframe(ChebyshevU(10),0.0..0.5))

    @test operator(dualdictionary(plat1,10,FourierWeight())) isa DiagonalOperator
    @test dualdictionary(plat2,10,measure(plat2)) isa ExtensionFrame
    @test operator(dualdictionary(plat3,10,ChebyshevWeight())) isa DiagonalOperator
    @test dualdictionary(plat4,10,measure(plat4)) isa ExtensionFrame
end
