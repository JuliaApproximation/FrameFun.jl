module test_suite_applications


using BasisFunctions
using FrameFuns
using FixedSizeArrays
using Base.Test
using Debug
FE = FrameFuns
BA = BasisFunctions

## Settings

# Test fourier extensions for all parameters
Extensive = false

# Show matrix vector product timings
const show_mv_times = true
total_mv_allocs = 0
total_mv_time = 0.0

const include_1d_tests = true
const include_2d_tests = true
const include_3d_tests = true
const include_bigfloat_tests = false

########
# Auxiliary functions
########

# Keep track of successes, failures and errors
global failures = 0
global successes = 0
global errors = 0

# Custom test handler
custom_handler(r::Test.Success) = begin print_with_color(:green, "#\tSuccess "); println("on $(r.expr)"); global successes+=1;  end
custom_handler(r::Test.Failure) = begin print_with_color(:red, "\"\tFailure "); println("on $(r.expr)\""); global failures+=1; end
custom_handler(r::Test.Error) = begin println("\"\t$(typeof(r.err)) in $(r.expr)\""); global errors+=1; end
#custom_handler(r::Test.Error) = Base.showerror(STDOUT,r);


# Check the accuracy of framefuns.
function msqerror_tol(f::Function, F; vals::Int=200, tol=1e-6)
    T = numtype(F)
    N = ndims(F)

    # Find the closest bounding grid around the domain
    TB = FE.boundingbox(FE.domain(F))

    point = Array{T}(N)
    elements = 0
    error = 0
    l = left(TB)
    r = right(TB)
    pvals_array = zeros(T,N)
    for i in 1:vals
        for i = 1:N
            pvals_array[i] = convert(T,rand())
        end
        pvals = Vec{N,T}(pvals_array)
        point = l+(r-l).*pvals
        if FE.in(point, FE.domain(F))
            elements += 1
            error += abs(f(point...)-F(point...))
            ## println("ratio ",f(point)/F(point...))
        end
    end
    @printf(" %3.2e",error/elements)
    return error>0 ? error/elements<tol : false
end

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

show_timings(F) = show_timings(F, F.approx_op)

show_timings(F, op::FE.FE_Solver) = show_timings(F, operator(op.problem))

# function show_timings(F, op::TensorProductOperator)
#     for i in 1:composite_length(op)
#         show_timings(F, element(op,i))
#     end
# end
#
# function show_timings(F, op::CompositeOperator)
#     for i in 1:composite_length(op)
#         show_timings(F, element(op,i))
#     end
# end
#
# show_timings(F, op::DimensionOperator) = show_timings(F, op.op)

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
    for Basis in (FourierBasis, ChebyshevBasis)
        println()
        println("## Basis: ", Basis)

        # Only 2 possible domains: an Interval and a Maskedgrid
        for D in [Interval(), Interval(-1.5,0.7), Interval(-1.5,-0.5)+Interval(0.5,1.5)]
            show(D); println()
            for n in [FE.default_frame_n(D, Basis) 99]
                println("\tN = $n")
                for n2 in (n-10, n, n+11)
                    # There is some symmetry around T=2, test smaller and larger values
                    for T in [1.7 FE.default_frame_T(D, Basis) 2.3]
                        print("T = $T \t")
                        B = Basis(n, -T, T)
                        F1 = Fun(f1, B, D)
                        F2 = Fun(f2, B, D)
                        F3 = Fun(f3, B, D) 
                        # Test Arithmetic
                        @test  msqerror_tol(x-> f1(x)+f2(x), F1+F2, tol=sqrt(eps(numtype(B)))*10)
                        @test  msqerror_tol(x-> f2(x)+f3(x), F2+F3, tol=sqrt(eps(numtype(B)))*10)
                        @test  msqerror_tol(x-> f1(x)*f2(x), F1*F2, tol=sqrt(eps(numtype(B)))*10)
                        @test  msqerror_tol(x-> f2(x)*f3(x), F2*F3, tol=sqrt(eps(numtype(B)))*10)
                        print("\t\t")
                    end
                end
            end
        end
        
    end
end

Test.with_handler(custom_handler) do
    test_arithmetics()
end

# Diagnostics
println()
println("Succes rate:\t$successes/$(successes+failures+errors)")
println("Failure rate:\t$failures/$(successes+failures+errors)")
println("Error rate:\t$errors/$(successes+failures+errors)")
if show_mv_times
    println("Total bytes in MV products:\t$total_mv_allocs")
    println("Total time in MV products:\t$total_mv_time")
end

(errors+failures)==0 || error("A total of $(failures+errors) tests failed")

end
