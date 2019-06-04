using FrameFun, DomainSets, Test, LinearAlgebra
@testset "Solvers" begin
    P = FourierExtensionPlatform(0.0..0.5)
    f = exp
    A,b = discretization(f, P, 100)
    M = FrameFun.plungeoperator(P,100)*A
    S = rSVD_solver(M;threshold=1e-4)
    x1 = S*FrameFun.plungeoperator(P,100)*b
    x2 = AZ_Zt(P,100)*(b-A*x1)
    @test norm(A*(x1+x2)-b) <  1e-3

    F = approximate(f, P, 100;REG=rSVD_solver)[1]
    x = .2; @test abs(F(x)-f(x)) < 1e-3

    P = FourierExtensionPlatform(0.0..0.5)
    f = exp
    A,b = discretization(f, P, 100)
    M = FrameFun.plungeoperator(P,100)*A
    S = pQR_solver(M;threshold=1e-4)
    x1 = S*FrameFun.plungeoperator(P,100)*b
    x2 = AZ_Zt(P,100)*(b-A*x1)
    @test norm(A*(x1+x2)-b) < 1e-3

    F = approximate(f, P, 100;REG=pQR_solver)[1]
    x = .2; @test abs(F(x)-f(x)) < 1e-3

    P = FourierExtensionPlatform(0.0..0.5)
    f = exp
    A,b = discretization(f, P, 100)
    M = FrameFun.plungeoperator(P,100)*A
    S = pSVD_solver(M;threshold=1e-4)
    x1 = S*FrameFun.plungeoperator(P,100)*b
    x2 = AZ_Zt(P,100)*(b-A*x1)
    @test norm(A*(x1+x2)-b) < 1e-3

    F = approximate(f, P, 100;REG=pSVD_solver)[1]
    x = .2; @test abs(F(x)-f(x)) < 1e-3
end
