module test_suite_applications


using BasisFunctions
using FrameFuns
using Base.Test
FE = FrameFuns
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
function test_arithmetics()
    f1(x) = 1
    f2(x) = sin(x^2-0.3)
    f3(x) = 1im*cos(x)
    @testset "Basis = $(name(instantiate(Basis,10))), D=$D" for Basis in (FourierBasis, ChebyshevBasis), D in [Interval(), Interval(-1.5,0.7), Interval(-1.5,-0.5)+Interval(0.5,1.5)]
        for n in [FE.default_frame_n(D, Basis) 99]
            for n2 in (n-10, n, n+11)
                # There is some symmetry around T=2, test smaller and larger values
                for T in [1.7 FE.default_frame_T(D, Basis) 2.3]
                    B = Basis(n, -T, T)
                    F1 = Fun(f1, B, D)
                    F2 = Fun(f2, B, D)
                    F3 = Fun(f3, B, D)
                    # Test Arithmetic
                    tol=sqrt(eps(numtype(B)))*10
                    @test  abserror(x-> f1(x)+f2(x), F1+F2) < tol
                    @test  abserror(x-> f2(x)+f3(x), F2+F3) < tol
                    @test  abserror(x-> f1(x)*f2(x), F1*F2) < tol
                    @test  abserror(x-> f2(x)*f3(x), F2*F3) < tol
                end
            end
        end
    end
end


test_arithmetics()


if show_mv_times
    println("Total bytes in MV products:\t$total_mv_allocs")
    println("Total time in MV products:\t$total_mv_time")
end

end
