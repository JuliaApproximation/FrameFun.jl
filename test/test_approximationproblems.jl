
using DomainSets, BasisFunctions, FrameFun
using Test

@testset "Platform approximation problems" begin
    ap1 = approximationproblem(Fourier(10))
    ap2 = approximationproblem(Fourier(10)^2)
    ap3 = approximationproblem(platform(Fourier(10)),10)
    ap4 = approximationproblem(Fourier(10), 0.0..0.5)
    ap5 = approximationproblem(Fourier(10), 20; samplingstyle=OversamplingStyle())

    @test length(factors(ap2))==2

    for ap in (ap1, ap2, ap3, ap4)
        @test !FrameFun.ap_hasproperty(ap, samplingparameter)
    end
    @test FrameFun.ap_hasproperty(ap5, samplingparameter)

    for ap in (ap1, ap2, ap3, ap4)
        if platform(ap) isa ProductPlatform
            FrameFun.cache!(ap, samplingparameter, (10,10))
        else
            FrameFun.cache!(ap, samplingparameter, 10)
        end
    end

    @test SamplingStyle(ap1) == InterpolationStyle()
    @test SamplingStyle(ap2) == ProductSamplingStyle(InterpolationStyle(),InterpolationStyle())
    @test SamplingStyle(ap3) == InterpolationStyle()
    @test SamplingStyle(ap4) == OversamplingStyle()
    @test SamplingStyle(ap5) == OversamplingStyle()

    @test SolverStyle(ap1) == InverseStyle()
    @test SolverStyle(ap2) == ProductSolverStyle(InverseStyle(),InverseStyle())
    @test SolverStyle(ap3) == InverseStyle()
    @test SolverStyle(ap4) == AZStyle()
    @test SolverStyle(ap5) == DirectStyle() # TODO: make this into an efficient least squares solver

    for ap in (ap1, ap2, ap3, ap4, ap5)
        @test FrameFun.ap_hasproperty(ap, samplingparameter)
    end

    @test dictionary(ap1) == Fourier(10)
    @test dictionary(ap2) == Fourier(10)^2
    @test dictionary(ap3) == Fourier(10)
    @test iscompatible(dictionary(ap4), extensionframe(Fourier(10), 0.0..0.5))
    @test dictionary(ap5) == Fourier(10)
end

@testset "Adaptive approximation problems" begin
    ap1 = approximationproblem(platform(Fourier(10)))
    ap2 = approximationproblem(platform(Fourier(10)^2))
    ap3 = approximationproblem(platform(extensionframe(Fourier(10), 0.0..0.5)))
    ap4 = approximationproblem(platform(Fourier(10)); samplingstyle=OversamplingStyle())

    for ap in (ap1,ap2,ap3,ap4)
        @test ap isa FrameFun.AdaptiveApproximation
    end

    @test dictionary(ap1,10)== Fourier(10)
    @test dictionary(ap2,10)== Fourier(10)^2
    @test iscompatible(dictionary(ap3,10),extensionframe(Fourier(10), 0.0..0.5))
    @test dictionary(ap4,10)== Fourier(10)

    @test SamplingStyle(ap1) == InterpolationStyle()
    @test SamplingStyle(ap4) == OversamplingStyle()
    ap1b = approximationproblem(ap1, 20)
    ap4b = approximationproblem(ap4, 20)
    @test ap1b isa FrameFun.PlatformApproximation
    @test ap4b isa FrameFun.PlatformApproximation
    @test SamplingStyle(ap1b) == InterpolationStyle()
    @test SamplingStyle(ap4b) == OversamplingStyle()
end
