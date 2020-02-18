using FrameFun, DomainSets, Test


@testset begin "Optimal adaptivity"
    P = FrameFun.ExtensionFramePlatform(WeightedSumPlatform(platform(ChebyshevT(10)^2), (x,y)->1.,
            (x,y)->sqrt(x^2+y^2)),.9*UnitDisk())
    f = (x,y) -> cos(pi*(x+y)) + sqrt(x^2+y^2)*sin(1+pi*(x+y))
    p = IncrementalCartesianParameterPath{2}()
    PP = parametrizedplatform(P, HierarchyPath(p,ProductPath(p,p)))

    F1, logbook1, _ = FrameFun.Adaptivity.adaptive_approximation(OptimalStyle(), f, PP; criterion = FNAStyle(),FNAerr=1e-2,FNAcoef=Inf)
    @test size(F1)==(45,)
    @test length(logbook1) == 8
    @test last(logbook1)[2] < .1
    @test norm(coefficients(F1))>10*2.188

    F, logbook, _ = FrameFun.Adaptivity.adaptive_approximation(OptimalStyle(), f, PP; criterion = FNAStyle(),FNAerr=1e-2,FNAcoef=10)
    @test size(F)==(190,)
    @test length(logbook) == 12
    @test last(logbook1)[2] < .1
    @test norm(coefficients(F))<10*2.188


    P = FrameFun.ExtensionFramePlatform(WeightedSumPlatform(parametrizedplatform(platform(ChebyshevT(10)^2),CartesianParameterPath((1,1))), (x,y)->1.,
            (x,y)->sqrt(x^2+y^2)),.9*UnitDisk())
    PP = parametrizedplatform(P)
    f = (x,y) -> cos(pi*(x+y)) + sqrt(x^2+y^2)*sin(1+pi*(x+y))


    F1 = Fun(f, PP; adaptivestyle = OptimalStyle(), criterion = FNAStyle(),FNAerr=1e-2,FNAcoef=Inf)
    @test size(F1)==(50,)
    F1 = Fun(f, PP; adaptivestyle = OptimalStyle(), criterion = FNAStyle(),FNAerr=1e-2,FNAcoef=10)
    @test size(F1)==(200,)
end
