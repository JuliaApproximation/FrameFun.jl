using FrameFun.ExtensionFrames
using BasisFunctions, Test, DomainSets
using BasisFunctions.Test: test_generic_dict_interface

import BasisFunctions.Test: point_outside_domain
using BasisFunctions.Test: generic_test_discrete_measure, generic_test_measure

point_outside_domain(::ExtensionFrame) = 1.

@testset "ExtensionFrame" begin
    d1 = extensionframe(Fourier(10),0.0..0.5)
    test_generic_dict_interface(d1)
    d2 = extensionframe(Fourier(10)^2,(0.0..0.5)^2)
    @test d1 == extensionframe(0.0..0.5, Fourier(10))
    @test d2 isa TensorProductDict
    @test d2 isa ExtensionFrameTensor
    @test basis(d1) == Fourier(10)
    @test basis(d2) == Fourier(10)^2
    @test support(d1) == 0.0..0.5
    @test support(d2) == (0.0..0.5)^2
    @test measure(d1) isa SubMeasure
    @test measure(d2) isa ProductMeasure
    @test iscompatible(extensiondual(d1, measure(d1)),element(extensiondual(d2, measure(d2)), 1))
end

@testset "SubMeasure" begin
    g = PeriodicEquispacedGrid(3,-1,1)^2
    d = UnitDisk()
    sg = subgrid(g,d)
    μ = discretemeasure(sg)
    @test μ isa DiscreteSubMeasure
    generic_test_discrete_measure(μ)
    @test supermeasure(μ) == discretemeasure(g)
    @test !(isnormalized(μ))
    @test restrict(supermeasure(μ),d)≈μ

    g = PeriodicEquispacedGrid(3,-1,1)×PeriodicEquispacedGrid(4,-1,1)
    d = UnitInterval()^2
    @assert subgrid(g,d) isa TensorSubGrid
    sg = subgrid(g,d)
    μ = discretemeasure(sg)
    generic_test_discrete_measure(μ )
    @test μ isa DiscreteTensorSubMeasure
    @test supermeasure(μ) ≈ discretemeasure(g)
    @test elements(μ) == (element(μ,1),element(μ,2))
    @test element(μ,1) == restrict(element(supermeasure(μ),1),UnitInterval())


    m = JacobiMeasure(rand(),rand())
    μ = restrict(m, UnitInterval())
    generic_test_measure(μ)
    @test supermeasure(μ) == m
    μ = submeasure(m, UnitInterval())
    generic_test_measure(μ)
    @test supermeasure(μ) == m
end
