
using DomainSets, BasisFunctions, FrameFun
using Test


@testset "ApproximationProblems" begin
    ap1 = approximationproblem(Fourier(10))
    ap2 = approximationproblem(Fourier(10)^2)
    ap3 = approximationproblem(platform(Fourier(10)),10)
    ap4 = approximationproblem(Fourier(10), 0.0..0.5)

    @test length(components(ap2))==2

    for ap in (ap1, ap2, ap3, ap4)
        @test samplingparam(ap) == nothing
    end

    for ap in (ap1, ap2, ap3, ap4)
        if ap isa ProductPlatformApproximation
            setsamplingparam!(ap, (10,10))
        else
            setsamplingparam!(ap, 10)
        end
    end

    @test SamplingStyle(ap1) == InterpolationStyle()
    @test SamplingStyle(ap2) == ProductSamplingStyle(InterpolationStyle(),InterpolationStyle())
    @test SamplingStyle(ap3) == InterpolationStyle()
    @test SamplingStyle(ap4) == OversamplingStyle()

    @test SolverStyle(SamplingStyle(ap1),ap1) == InverseStyle()
    @test SolverStyle(SamplingStyle(ap2),ap2) == ProductSolverStyle(InverseStyle(),InverseStyle())
    @test SolverStyle(SamplingStyle(ap3),ap3) == InverseStyle()
    @test SolverStyle(SamplingStyle(ap4),ap4) == DirectStyle()

    for ap in (ap1, ap2, ap3, ap4)
        @test samplingparam(ap) != nothing
    end

    @test dictionary(ap1) == Fourier(10)
    @test dictionary(ap2) == Fourier(10)^2
    @test dictionary(ap3) == Fourier(10)
    @test iscompatible(dictionary(ap4), extensionframe(Fourier(10), 0.0..0.5))

    ap1 = approximationproblem(platform(Fourier(10)))
    ap2 = approximationproblem(platform(Fourier(10)^2))
    ap3 = approximationproblem(platform(extensionframe(Fourier(10), 0.0..0.5)))

    for ap in (ap1,ap2,ap3)
        @test ap isa AdaptiveApproximation
    end

    @test dictionary(ap1,10)== Fourier(10)
    @test dictionary(ap2,10)== Fourier(10)^2
    @test iscompatible(dictionary(ap3,10),extensionframe(Fourier(10), 0.0..0.5))
end
