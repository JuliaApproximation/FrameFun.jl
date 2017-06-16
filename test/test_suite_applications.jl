module test_suite_applications

using Domains
using BasisFunctions
using FrameFun
using Base.Test
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

show_timings(F) = show_timings(F, F.approx_op)

show_timings(F, op::FE.FE_Solver) = show_timings(F, operator(op.problem))

function show_timings(F, op)
    if show_mv_times
        c1 = zeros(eltype(op), size(src(op)))
        c2 = zeros(eltype(op), size(dest(op)))
        apply!(op, c2, c1)
        t = @timed apply!(op, c2, c1)
        if t[3] > 700
            print_with_color(:red, "\"\tMV product: ")
        elseif t[3] == 0
            print_with_color(:green, "\"\tMV product: ")
        else
            print_with_color(:blue, "\"\tMV product: ")
        end
        println(t[3], " bytes allocated, in ", t[2], " seconds")
        # end
        global total_mv_allocs += t[3]
        global total_mv_time += t[2]
    end
end

#######
# Testing -- Accuracy
#######



function test_differential_equations_1d()
    @testset "diff 1D" for solver in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
        B = FourierBasis(101,-1,1)
        Dom = interval(-0.5,0.5)
        # Set up Boundary conditions
        diff = IdentityOperator(B)
        df(x) = 0;
        BC = BoundaryCondition(B,diff,Dom,df)
        # Set up Differential equation
        f(x) = x
        Diff = differentiation_operator(B)^2
        DE = DiffEquation(B,Dom,Diff, f, (BC,))
        # Actually solve the differential equation
        F = solve(DE, solver=solver)
        sol(x) = x^3/6 - x/24
        error = abserror(sol,F)
        @test (error < sqrt(eps(numtype(B)))*10)
        error = abserror(f,F'')
        @test (error < sqrt(eps(numtype(B)))*100)
    end
end

function test_differential_equations_2d()
    @testset "diff 2D" for solver in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
        B = FourierBasis(11,-1,1)⊗FourierBasis(11,-1,1)
        Dom = disk(0.8)
        # Set up Boundary conditions
        diff = IdentityOperator(B)
        df(x,y) = x-y;
        BC = BoundaryCondition(B,diff,Dom,df)
        # Set up Differential equation
        f(x,y) = 0
        Diff = differentiation_operator(B,(2,0))+differentiation_operator(B,(0,2))
        DE = DiffEquation(B,Dom,Diff, f, (BC,))
        # Actually solve the differential equation
        F = solve(DE, solver=solver)
        error = abserror(df,F)
        @test (error < 0.3)
        error = abserror(f,∂x(∂x(F))+∂y(∂y(F)))
        @test (error < 0.3)
    end
end

# We just test the accuracy, not the smoothing properties
function test_smoothing_1d()
    @testset "Smoothing $(name(instantiate(Basis,10)))" for Basis in (FourierBasis, ChebyshevBasis)
        B = Basis(101,-1,1)
        D = interval(-0.5,0.5)
        f(x) = exp(x)
        fscale(i) = 10.0^-4+abs(i)+abs(i)^2+abs(i)^3
        F = Fun(f,B,D;solver=FrameFun.FE_SmoothProjectionSolver,scale=fscale)
        # F = Fun(f,B,D;solver=FrameFun.FE_ProjectionSolver)
        @test (abserror(f,F) < 100*sqrt(eps(numtype(B))))
    end
end

function test_smoothing_2d()
    @testset "Smoothing $(name(instantiate(Basis,10)))" for Basis in (FourierBasis, ChebyshevBasis)
        B = Basis(20,-1.0,1.0)⊗Basis(20,-1.0,1.0)
        D = disk(0.5)
        f(x,y) = exp(x*y)
        fscale(i,j) = 10.0^-4+100*abs((i)^2+abs(j^2))
        F = Fun(f,B,D;solver=FrameFun.FE_SmoothProjectionSolver,scale=fscale)
        @test (abserror(f,F) < 100*sqrt(sqrt(eps(numtype(B)))))
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
