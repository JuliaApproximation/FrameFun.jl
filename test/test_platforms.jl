using Test, FrameFun.Platforms, BasisFunctions

@testset "Generic Platforms" begin
    P = platform(Fourier(10))
    @test dictionary(P,param(Fourier(10))) == Fourier(10)
    @test P[10] == Fourier(10)
    @test operator(dualdictionary(P,10,measure(P))) ≈ operator(dual(Fourier(10)))
    @test DictionaryStyle(P)==BasisStyle()
    @test SamplingStyle(P) == InterpolationStyle()
    @test SolverStyle(P,SamplingStyle(P)) == InverseStyle()
    @test param_first(P) == 10
    @test param_double(P,10) == 20
    @test param_increment(P,10) == 11
    @test param_inbetween(P,10,20) == 15
    @test correctparamformat(P,10)
    @test !correctparamformat(P,(10,))

    @test elements(platform(Fourier(10)^2)) == (P,P)
    @test element(platform(Fourier(10)^2),1) == P

    P = platform(Fourier(10)^2)
    @test correctparamformat(P,(10,10))
    @test !correctparamformat(P,(10,))
    @test !correctparamformat(P,10)
    @test param_first(P) == (10,10)
    @test param_double(P,(10,20)) == (20,40)
    @test param_increment(P,(10,11)) == (11,12)
    @test param_inbetween(P,(10,10),(20,10)) == (15,10)
    @test ProductPlatform(platform(Fourier(10)),2) == P
    @test element(P,1) == platform(Fourier(10))
    length(elements(P))==2
    @test dictionary(P,10) == Fourier(10)^2
    @test dictionary(P,param(Fourier(10)^2)) == Fourier(10)^2
    @test P[10] == Fourier(10)^2
    @test dictionary(P,(11,12)) == Fourier(11)⊗Fourier(12)
    @test dualdictionary(P,(10,10), measure(Fourier(10)^2)) isa TensorProductDict
    isbasis(Fourier(10)^2)
    @test DictionaryStyle(P)==BasisStyle()
    @test SamplingStyle(P) == ProductSamplingStyle(InterpolationStyle(),InterpolationStyle())
    @test SolverStyle(P,SamplingStyle(P)) == ProductSolverStyle(InverseStyle(),InverseStyle())
end

@testset "Platform Styles" begin
    @test SamplingStyle(platform(Fourier(10))) == InterpolationStyle()
    @test SolverStyle(platform(Fourier(10)), InterpolationStyle()) == InverseStyle()
    @test SolverStyle(platform(Fourier(10)), OversamplingStyle()) == DirectStyle()
    @test SolverStyle(platform(Fourier(10)), GridStyle()) == DirectStyle()
    @test SolverStyle(platform(Fourier(10)), GenericSamplingStyle()) == DirectStyle()


    @test SolverStyle(FrameStyle(), platform(Fourier(10)), InterpolationStyle()) == AZStyle()
    @test SolverStyle(FrameStyle(), platform(Fourier(10)), OversamplingStyle()) == AZStyle()
    @test SolverStyle(FrameStyle(), platform(Fourier(10)), GridStyle()) == AZStyle()
    @test SolverStyle(FrameStyle(), platform(Fourier(10)), GenericSamplingStyle()) == AZStyle()

    @test DictionaryStyle(platform(Fourier(10))) == BasisStyle()
    @test DictionaryStyle(platform(MultiDict((Fourier(10),Fourier(10))))) == UnknownDictionaryStyle()
    @test elements(SolverStyle(platform(Fourier(10)^2), SamplingStyle(platform(Fourier(10)^2)))) == (InverseStyle(), InverseStyle())
    @test ProblemStyle(platform(Fourier(10)))  == DictionaryOperatorStyle()

end
