using FrameFun.ExtensionFrames
using BasisFunctions, Test, DomainSets
using BasisFunctions.Test: test_generic_dict_interface

import BasisFunctions.Test: point_outside_domain

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
