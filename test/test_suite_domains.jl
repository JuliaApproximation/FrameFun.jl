using DomainSets, BasisFunctions, FrameFun

using Test

FE = FrameFun
BA = BasisFunctions

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

function test_extended_domains_3d()
    domains = (atomium(),)
end

function test_extended_domains_2d()
    D1 = characteristic(x->(x[1]-x[2])<0,disk(1.0))
    @test in([0.5;1.0],D1)
    @test ~in([1.0,0.5],D1)
    basis = ChebyshevBasis(21)⊗ChebyshevBasis(21)
    domain = disk(1.0)
    f = (x,y)->exp(x+y)
    F = Fun(f, basis, domain)
    F2 = F*2
    D2 = F<2
    D3 = F<F2
    @test ~in([0.5;0.5],D2)
    @test in([0.25;0.25],D2)
    @test in([0.5;0.5],D3)
    D4 = mandelbrot()
    @test in([0.0;0.0],D4)
    @test ~in([-0.5;0.5],D4)
    D5 = juliaset()
    @test in([0.6;0.0],D5)
    @test ~in([-0.1;-0.3],D5)
    D6 = polardomain(x->1+0.2*cos(5x),disk(2.0))
    @test ~in([1.25;0],D6)
    @test in([0.79;0],D6)
end

domn(point,domain::Domain,normalvector; tol=1e-14) = norm(normal(point,domain)-normalvector)<tol
domb(point,domain::Domain; tol=1e-14) = abs(distance(point,domain))<tol

function test_distances_and_normals()
    D1 = simplex(Val{2})
    @test domb([0.5,0.5],D1)
    @test domn([0.5,0.5],D1,[1/sqrt(2),1/sqrt(2)])
    @test domb([0.5,0],D1)
    @test domn([0.5,0.0],D1,[0,-1])

    D6 = polardomain(x->1+0.2*cos(5x),disk(2.0))
    @test domb([-0.8,0],D6)
    @test domn([-0.8,0],D6,[-1,0],tol=1e-12)

    D2 = disk(sqrt(0.5))
    @test domb([0.5,0.5],D2)
    @test domn([0.5,0.5],D2,[1/sqrt(2),1/sqrt(2)])

    D2b = circle(sqrt(0.5))
    @test domb([0.5,0.5],D2b)
    @test domn([0.5,0.5],D2b,[1/sqrt(2),1/sqrt(2)])

    D3 = -1..1.0
    @test isapprox(distance(0.5,D3),0.5)
    @test isapprox(normal(1.0,D3),1)

    D5 = cube(Val{2})
    @test domb([0.4,0.0],D5)
    @test domn([0.4,0.0],D5,[0,-1])

    D4 = UnionDomain(D6, disk())\D2
    @test domb([sqrt(0.5),0.0],D4)
    @test domb([-1.0,0.0],D4)
    @test ~domb([-0.8,0.0],D4)
    @test domn([-1.0,0],D4,[-1,0])
    @test domn([0.5,0.0],D4,[-1.0,0])

end



delimit("Extended Domains")

test_extended_domains_3d()

test_extended_domains_2d()

test_distances_and_normals()
