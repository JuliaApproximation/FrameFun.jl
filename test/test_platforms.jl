module test_suite_applications

using BasisFunctions, DomainSets, FrameFun
using Test


@testset "WeightedSumPlatform" begin
    # Simple platform construction test
    n = 12
    domain = 0.0..0.5
    P = FourierExtensionPlatform(domain)
    WP = WeightedSumPlatform(P,[x->sqrt(x),x->1])
    f = x->sqrt(x)*(1-x)-exp(x)
    F = Fun(f, WP, n, solverstyle=AZStyle())
    rgrid = randomgrid(domain, 200)
    abserror = sum(abs.(F.(rgrid)-f.(rgrid)))/length(rgrid)
    @test abserror<1e-7
end

end # module
