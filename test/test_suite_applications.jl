module test_suite_applications

using DomainSets, BasisFunctions, FrameFun
using Test

FE = FrameFun
BA = BasisFunctions

## Settings

# Test fourier extensions for all parameters
Extensive = false

# Show matrix vector product timings
const show_mv_times = false
total_mv_allocs = 0
total_mv_time = 0.0

########
# Auxiliary functions
########

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end


#######
# Testing -- Accuracy
#######



function test_differential_equations_1d()
    @testset "diff 1D dirichlet" for solverstyle in (AZStyle(), DirectStyle())
        B = Fourier(101) → -1..1
        Dom = Interval(-0.5,0.5)
        # Set up Boundary conditions
        BC = DirichletBC(x->0)
        # Set up Differential equation
        f(x) = x
        Diff = differentiation(B)*differentiation(B)
        DE = DiffEquation(B, Dom, Diff, f, (BC,))
        # Actually solve the differential equation
        F = solve(DE, solverstyle=solverstyle, directsolver=:qr)
        sol(x) = x^3/6 - x/24
        error = abserror(sol,F)
        @test (error < sqrt(eps(real(codomaintype(B))))*10)
        error = abserror(f,F'')
        @test (error < sqrt(eps(real(codomaintype(B))))*100)
    end
    @testset "diff 1D mixed" for solverstyle in (AZStyle(), DirectStyle())
        B = Fourier(101) → -2..2
        Dom = Interval(-1.0,1.0)
        # Set up Boundary conditions
        BC1 = DirichletBC(x->0, Interval(-2.0,0.0))
        BC2 = NeumannBC(x->0, Interval(-2.0,0.0))
        # Set up Differential equation
        f(x) = x
        Diff = differentiation(B)*differentiation(B)
        DE = DiffEquation(B, Dom, Diff, f, (BC1,BC2))
        # Actually solve the differential equation
        F = solve(DE, solverstyle=solverstyle, directsolver=:qr)
        sol(x) = x^3/6 - x/2 -1/3
        error = abserror(sol,F)
        @test (error < sqrt(eps(real(codomaintype(B))))*10)
        error = abserror(f,F'')
        @test (error < sqrt(eps(real(codomaintype(B))))*100)
    end
end

function test_differential_equations_2d()
    # @testset "diff 2D dirichlet" for solverstyle in (AZStyle(), DirectStyle())
    @testset "diff 2D dirichlet" for solverstyle in (AZStyle(),)
        B = (Fourier(11) → -1..1)^2
        Dom = disk(0.8)
        # Set up Boundary conditions
        df = (x,y)->x-y
        BC = DirichletBC(df,DomainSets.euclideanspace(Val{2}()))
        # Set up Differential equation
        f(x,y) = 0
        Diff = differentiation(B,(2,0))+differentiation(B,(0,2))
        DE = DiffEquation(B, Dom, Diff, f, (BC,))
        # Actually solve the differential equation
        F = solve(DE, solverstyle=solverstyle, directsolver=:qr)
        error = abserror(df,F)
        @test (error < 0.3)
        error = abserror(f,∂x(∂x(F))+∂y(∂y(F)))
        @test (error < 0.3)
    end
    # @testset "diff 2D Neumann" for solverstyle in (AZStyle(), DirectStyle())
    # @testset "diff 2D Neumann" for solverstyle in (AZStyle(),)
    #     B = (Fourier(11) → -1..1)^2
    #     Dom = disk(0.8)
    #     # Set up Boundary conditions
    #     df = (x,y)->x-y
    #     BC = NeumannBC(df,DomainSets.euclideanspace(Val{2}()))
    #     # Set up Differential equation
    #     f(x,y) = 0
    #     Diff = differentiation(B,(2,0))+differentiation(B,(0,2))
    #     DE = DiffEquation(B, Dom, Diff, f, (BC,))
    #     # Actually solve the differential equation
    #     F = solve(DE, solverstyle=solverstyle, directsolver=:qr)
    #     # We should find a way to check the boundary condition more easily
    #     berror = sum(abs.(element(operator(DE),2,1)*coefficients(Diff*F)-element(FrameFun.rhs(DE),2)))/length(element(FrameFun.rhs(DE),2))
    #     @test berror < 0.6
    #     error = abserror(f,∂x(∂x(F))+∂y(∂y(F)))
    #     @test (error < 0.3)
    # end
end

# We just test the accuracy, not the smoothing properties
function test_smoothing_1d()
    @testset "Smoothing $(name(Basis(10)))" for Basis in (ChebyshevT,)
        basis = Basis(101,-1,1)
        dom = -0.5..0.5
        f = exp
        fscale(dict, i) = 10.0^-4+i+i^2+i^3
        F = Fun(f, basis, dom; solverstyle = AZSmoothStyle(), scaling=fscale)
        @test (abserror(f,F) < 100*sqrt(eps(Float64)))
    end
end

function test_smoothing_2d()
    @testset "Smoothing $(name(Basis(10)))" for Basis in (ChebyshevT,)
        basis = Basis(20,-1.0,1.0)⊗Basis(20,-1.0,1.0)
        dom = disk(0.5)
        f = (x,y) -> exp(x*y)
        fscale(dict, I) = 10.0^-4+100*(I[1]^2+I[2]^2)
        F = Fun(f, basis, dom; solverstyle = AZSmoothStyle(), scaling=fscale)
        @test (abserror(f,F) < 100*sqrt(sqrt(eps(Float64))))
    end
end


test_smoothing_1d()

test_smoothing_2d()

test_differential_equations_1d()

test_differential_equations_2d()




if show_mv_times
    println("Total bytes in MV products:\t$total_mv_allocs")
    println("Total time in MV products:\t$total_mv_time")
end

end
