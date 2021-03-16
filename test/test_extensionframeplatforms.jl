using BasisFunctions, FrameFun.ExtensionFrames, FrameFun.ExtensionFramePlatforms,
    DomainSets, Test, StaticArrays, FrameFun.Platforms, FrameFun.FrameFunInterface


@testset "ExtensionFramePlatforms" begin
    d1 = 0.0..0.5
    d2 = (0.0..0.5)^2
    d3 = disk(.4,SVector(.5,.5))
    dict1 = extensionframe(Fourier(10),d1)
    dict2 = extensionframe(Fourier(10)^2,d2)
    dict3 = extensionframe(Fourier(10)^2,d3)
    plat1 = (platform(dict1),10)
    plat2 = (platform(dict2),(10,10))
    plat3 = (platform(dict3),(10,10))

    @test dictionary(plat1...) isa ExtensionFrame
    @test dictionary(plat2...) isa ExtensionFrameTensor
    @test dictionary(plat3...) isa ExtensionFrame

    @test oversampling_grid(plat1...;oversamplingfactor=2)≈subgrid(FourierGrid(38),d1)
    @test oversampling_grid(plat2...;oversamplingfactor=2)≈subgrid(FourierGrid(27)^2,d2)
    @test oversampling_grid(plat3...;oversamplingfactor=2)≈subgrid(FourierGrid(21)^2,d3)
    @test length(oversampling_grid(plat1...;oversamplingfactor=2)) == 20
    @test abs(length(oversampling_grid(plat2...;oversamplingfactor=2)) - 200)<20
    @test abs(length(oversampling_grid(plat3...;oversamplingfactor=2)) - 200 )< 20

    @test SamplingStyle(platform(dict1)) == OversamplingStyle()
    @test SamplingStyle(platform(dict2)) == ProductSamplingStyle(OversamplingStyle(),OversamplingStyle())
    @test SamplingStyle(platform(dict3)) == OversamplingStyle()

    @test SolverStyle(platform(dict1),OversamplingStyle()) == AZStyle()
    @test SolverStyle(platform(dict2),ProductSamplingStyle(OversamplingStyle(),OversamplingStyle())) == ProductSolverStyle(AZStyle(),AZStyle())
    @test SolverStyle(platform(dict3),OversamplingStyle()) == AZStyle()

    @test DictionaryStyle(platform(dict1)) == FrameStyle()
    @test DictionaryStyle(platform(dict2)) == FrameStyle()
    @test DictionaryStyle(platform(dict3)) == FrameStyle()

    ddict1 = dualdictionary(plat1...,measure(plat1[1]))
    @test ddict1 isa ExtensionFrame
    @test operator(basis(ddict1)) ≈ operator(dual(Fourier(10),FourierWeight()))
    ddict2 = dualdictionary(plat2...,measure(plat2[1]))
    @test ddict2 isa ExtensionFrameTensor
    ddict3 = dualdictionary(plat3...,measure(plat3[1]))
    @test ddict3 isa ExtensionFrame
end
