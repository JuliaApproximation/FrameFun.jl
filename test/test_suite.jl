module test_suite


using BasisFunctions
using FrameFun
using StaticArrays
using Base.Test
FE = FrameFun
BA = BasisFunctions

## Settings

# Test fourier extensions for all parameters
Extensive = false

# Show matrix vector product timings
const show_mv_times = false
const verbose = true
total_mv_allocs = 0
total_mv_time = 0.0

const include_1d_tests = true
const include_2d_tests = true
const include_3d_tests = true
const include_bigfloat_tests = false

########
# Auxiliary functions
########



function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

show_timings(F) = show_mv_times && show_timings(F, F.approx_op)

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

function test_1d_cases()
    delimit("1D")

    # Complex and real functions / Float64
    f(x) = cos(x.^2-0.5)-1
    g(x) = 1im*cos(x.^2-0.5)-1
    # Chebyshev and Fourier Bases

    @testset "result" for ELT in (Float32,Float64), Basis in (FourierBasis, ChebyshevBasis), D in [Interval(), Interval(-1.5,0.7), Interval(-1.5,-0.5)+Interval(0.5,1.5)], solver in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
        println()
        println("Testing \t solver = $solver, \n\t\t Domain = $D, \n\t\t Basis = $(name(instantiate(Basis,10))),\n\t\t ELT = $ELT ")
        verbose && println("N\t T\t Complex?\t abserror\t time\t\t \memory   ")

        for n in [FE.default_frame_n(D, Basis) 99]

            # There is some symmetry around T=2, test smaller and larger values
            for T in [1.7 FE.default_frame_T(D, Basis) 2.3]
                for func in (f,g)

                    B = Basis(n, -T, T, ELT)
                    F = @timed( Fun(func, B, D; solver=solver) )
                    error = abserror(func, F[1])
                    if verbose
                        print("$n\t $T\t\t")
                        func==g ? print("Y\t") : print("N\t")
                        @printf("%3.2e\t %3.2e s\t %3.2e bytes \n",error, F[2],F[3])
                    end
                    @test  (error < sqrt(eps(ELT))*10)
                    show_timings(F[1])
                end
            end
        end
    end
end

function test_bigfloat()
    f(x) = cos(x.^2) - big(1.0)
    g(x) = big(1.0)im * cos(x.^2) - big(1.0)
    @testset "result" for Basis in (FourierBasis, ChebyshevBasis),  D in [Interval(BigFloat(-3//2),BigFloat(7//10))]
        println()
        println("Testing \t solver = FE.FE_DirectSolver{ELT}\n\t\t Domain = $D, \n\t\t Basis = $(name(instantiate(Basis,10))),\n\t\t ELT = BigFloat ")
        verbose && println("N\t T\t Complex?\t abserror\t time\t\t \memory   ")
        for T in [BigFloat(17//10) FE.default_frame_T(D, Basis) BigFloat(23//10)]
            for func in (f,g)
                B = Basis(91, -T, T)
                F = @timed( Fun(func, B, D; solver = FE.FE_DirectSolver) )
                error = abserror(func, F[1])
                if verbose
                    @printf("91\t %3.2e\t\t",T)
                    func==g ? print("Y\t") : print("N\t")
                    @printf("%3.2e\t %3.2e s\t %3.2e bytes \n",error, F[2],F[3])
                end
                @test  error < 1e-20
                show_timings(F)
            end
        end
    end
end

function test_2d_cases()
    delimit("2D")

    f(x,y) = cos(0.5*x)+2*sin(0.2*y)-1.0*x*y
    g(x,y) = 1im*cos(0.5*x)+2*sin(0.2*y)-1.0im*x*y
    @testset "result" for Basis in (FourierBasis, ChebyshevBasis), D in [Disk(), Disk(1.2,[-0.1,-0.2]), Cube((-1.0,-1.5),(0.5,0.7))], solver in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
        println()
        println("Testing \t solver = $solver \n\t\t Domain = $D, \n\t\t Basis = $(name(instantiate(Basis,10)⊗instantiate(Basis,10))),\n\t\t ELT = Float64 ")
        verbose && println("N\t\t T\t\t Complex?\t abserror\t time\t\t \memory   ")

        for n in ((11,11),)
            for T in (Extensive ? (FE.default_frame_T(D, Basis),) : ((1.7,1.7),FE.default_frame_T(D, Basis),(2.3,2.3)))

                B = Basis(n[1],-T[1],T[1]) ⊗ Basis(n[2],-T[2],T[2])
                for func in (f,g)
                    F = @timed( Fun(func, B, D; solver=solver))
                    error = abserror(func, F[1])
                    if verbose
                        print("$n \t $T \t\t")
                        func==g ? print("Y\t") : print("N\t")
                        @printf("%3.2e\t %3.2e s\t %3.2e bytes \n",error, F[2],F[3])
                    end
                    @test  error < 1e-3
                    show_timings(F)
                end

            end
        end
    end
end

function test_3d_cases()
    delimit("3D")

    f(x,y,z) = cos(x)+sin(y)-x*z
    @testset "result" for Basis in (FourierBasis, ChebyshevBasis), D in (Cube((-1.2,-1.0,-0.9),(1.0,0.9,1.2)),FE.tensorproduct(Interval(-1.0,1.0),Disk(1.05)), FE.Ball(1.2,[-0.3,0.25,0.1])), solver in (FE.FE_ProjectionSolver, )
                show(solver); println()
        println()
        println("Testing \t solver = $solver \n\t\t Domain = $D, \n\t\t Basis = $(name(instantiate(Basis,10)⊗instantiate(Basis,10))),\n\t\t ELT = Float64 ")
        verbose && println("N\t\t T\t\t Complex?\t abserror\t time\t\t \memory   ")

        n = FE.default_frame_n(D, Basis)

        for T in ((2.3,2.4,1.7), FE.default_frame_T(D, Basis))
            B = Basis(n[1],-T[1],T[1]) ⊗ Basis(n[2],-T[2],T[2]) ⊗ Basis(n[3],-T[3],T[3])
            F = @timed( Fun(f, B, D; solver=solver, cutoff=10.0^(3/4*log10(eps(numtype(B)))),sampling_factor=1.5))
            error = FrameFun.residual(f, F[1])/length(B)
            if verbose
                print("$n \t $T \t\t")
                print("N\t")
                @printf("%3.2e\t %3.2e s\t %3.2e bytes \n",error, F[2],F[3])
            end
            @test  error < 1e-2
            show_timings(F)
        end
    end
end


    delimit("Algorithm Implementation and Accuracy")
    # FFTW.set_num_threads(Sys.CPU_CORES)

    if include_1d_tests
        test_1d_cases()
    end

    if include_bigfloat_tests
        test_bigfloat()
    end

    if include_2d_tests
        test_2d_cases()
    end

    if include_3d_tests
        test_3d_cases()
    end

    delimit("Random circles")
    dom = FE.randomcircles(10)
#    b = FourierBasis(21) ⊗ FourierBasis(21)
    b = FourierBasis(21) ⊗ ChebyshevBasis(21)
    f(x,y) = cos(20*x+22*y)
    @time F = Fun(f,b,dom)
    show_timings(F)

if show_mv_times
    println("Total bytes in MV products:\t$total_mv_allocs")
    println("Total time in MV products:\t$total_mv_time")
end


end
