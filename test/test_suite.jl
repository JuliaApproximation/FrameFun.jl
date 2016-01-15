module test_suite


using BasisFunctions
using FrameFuns
using FixedSizeArrays
using Base.Test
FE = FrameFuns
BA = BasisFunctions

## Settings

# Test fourier extensions for all parameters
Extensive = false
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
    N = dim(F)

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
    println("############")
    println("# ",s)
    println("############")
end

#######
# Testing
#######


Test.with_handler(custom_handler) do


    delimit("Algorithm Implementation and Accuracy")
    delimit("1D")

    # Complex and real functions / Float64
    f(x) = cos(x.^2)-1.0
    g(x) = 1im*cos(x.^2)-1.0
    # Chebyshev and Fourier Bases
    for Basis in (FourierBasis, ChebyshevBasis)
        println()
        println("## Basis: ", Basis)

        # Only 2 possible domains: an Interval and a Maskedgrid
        for D in [Interval(-1.5,0.7) Interval(-1.5,-0.5)+Interval(0.5,1.5)]
            show(D); println()

            # 2 possible solvers
            for solver in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
                show(solver); println()

                for n in [FE.default_frame_n(D, Basis) 201]
                    println("\tN = $n")

                    # There is some symmetry around T=2, test smaller and larger values
                    for T in [1.7 FE.default_frame_T(D, Basis) 2.3]
                        print("T = $T \t")

                        for func in (f,g)
                            F = @timed( Fun(Basis, func, D; solver=solver, n=n, T=T) )
                            F = @timed( Fun(Basis, func, D; solver=solver, n=n, T=T) )

                            @printf("%3.2e s\t %3.2e bytes",F[2],F[3])
                            @test  msqerror_tol(func, F[1], tol=1e-7)
                            if func == f
                                print("\t\t")
                            end
                        end
                    end
                end
            end
        end
    end

    ## Float32
    f(x) = cos(x.^2)-1.0f0
    g(x) = 1.0f0im*cos(x.^2)-1.0f0
    # Chebyshev and Fourier Bases
    for Basis in (FourierBasis, ChebyshevBasis)
        println()
        println("## Basis: ", Basis)
        # Only 2 possible domains: an Interval and a Maskedgrid
        for D in [Interval(-1.5f0,0.7f0) Interval(-1.5f0,-0.5f0)+Interval(0.5f0,1.5f0)]
            show(D); println()

            # 2 possible solvers
            for solver in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
                show(solver); println()

                for n in [FE.default_frame_n(D, Basis) 99]
                    println("\tN = $n")
                    # There is some symmetry around T=2, test smaller and larger values
                    for T in [1.7 FE.default_frame_T(D, Basis) 2.3]
                        print("T = $T \t")

                        for func in (f,g)

                            F = @timed( Fun(Basis, func, D; solver=solver, n=n, T=T))
                            F = @timed( Fun(Basis, func, D; solver=solver, n=n, T=T))

                            @printf("%3.2e s\t %3.2e bytes",F[2],F[3])
                            @test  msqerror_tol(func,F[1],tol=1e-4)
                            if func==f
                                print("\t\t")
                            end
                        end
                    end
                end
            end
        end
    end

    test_bigfloat = true
    if test_bigfloat
        f(x) = cos(x.^2) - big(1.0)
        g(x) = big(1.0)im * cos(x.^2) - big(1.0)
        for Basis in (FourierBasis, ChebyshevBasis)
            println()
            println("## Basis: ", Basis)
            # Some BigFloat tests (DirectSolver only)
            for D in [Interval(BigFloat(-3//2),BigFloat(7//10)) Interval(BigFloat(-3//2),BigFloat(-1//2))+Interval(BigFloat(1//2),BigFloat(3//2))]
                show(D); print("\n")
                for T in [BigFloat(17//10) FE.default_frame_T(D, Basis) BigFloat(23//10)]
                    print("T = $T \t")
                    for func in (f,g)

                        F = @timed( Fun(Basis, func, D; solver = FE.FE_DirectSolver, n=71, T=T) )

                        @printf("%3.2e s\t %3.2e bytes",F[2],F[3])
                        @test  msqerror_tol(func,F[1],tol=1e-20)
                        if func==f
                            print("\t\t")
                        end
                    end
                end
            end
        end
    end

    delimit("2D") 

    f(x,y) = cos(0.5*x)+2*sin(0.2*y)-1.0*x*y
    g(x,y) = 1im*cos(0.5*x)+2*sin(0.2*y)-1.0im*x*y
    for Basis in (FourierBasis, ChebyshevBasis)
        println()
        println("## Basis: ", Basis)

        for D in [Disk(2.0,[-2.0,-2.0]) Cube((-1.0,-1.5),(0.5,0.7))]
            show(D); println()

            for solver in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
                show(solver); println()

                for n in ((11,11),)
                    println("\tN = $n")
                    for T in (Extensive ? (FE.default_frame_T(D, Basis),) : ((1.7,1.7),FE.default_frame_T(D, Basis),(2.3,2.3)))
                        print("T = $T\t")
                        for func in (f,g)

                            F = @timed( Fun(Basis, func, D; solver=solver, n=n, T=T))

                            @printf("%3.2e s\t %3.2e bytes",F[2],F[3])
                            @test msqerror_tol(func, F[1], tol=1e-3)
                            if func==f
                                print("\t\t")
                            end
                        end
                    end
                end
            end
        end
    end

    delimit("3D")

    f(x,y,z) = cos(x)+sin(y)-x*z
    for Basis in (FourierBasis, ChebyshevBasis)
        println()
        println("## Basis: ", Basis)

        for D in (Cube((-0.2,-2.0,0.0),(1.0,-0.1,1.2)),FE.TensorProductDomain(Interval(-1.0,1.0),Disk(1.05)), FE.Ball(1.2,[-1.3,0.25,1.0]))
            show(D); println()
            for solver in (FE.FE_ProjectionSolver, )
                show(solver); println()

                n = FE.default_frame_n(D, Basis)
                println("\tN = $n")

                for T in ((1.7,1.7,1.7), FE.default_frame_T(D, Basis))
                    print("T = $T\t")
                    F = @timed( Fun(Basis, f, D; solver=solver, n=n, T=T))
                    @printf("%3.2e s\t %3.2e bytes", F[2], F[3])
                    @test msqerror_tol(f, F[1], tol=1e-2)
                end
            end
        end
    end
end

# Diagnostics
println()
println("Succes rate:\t$successes/$(successes+failures+errors)")
println("Failure rate:\t$failures/$(successes+failures+errors)")
println("Error rate:\t$errors/$(successes+failures+errors)")

(errors+failures)==0 || error("A total of $(failures+errors) tests failed")

end

