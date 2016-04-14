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
    println()
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
    f(x) = cos(x.^2-0.5)-1
    g(x) = 1im*cos(x.^2-0.5)-1
    # Chebyshev and Fourier Bases
    
    for ELT in (Float32,Float64)
        for Basis in (FourierBasis, ChebyshevBasis)
            println()
            println("## Basis: ", Basis)

            # Only 2 possible domains: an Interval and a Maskedgrid
            for D in [Interval(), Interval(-1.5,0.7), Interval(-1.5,-0.5)+Interval(0.5,1.5)]
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
                                B = Basis(n, -T, T, ELT)
                                F = @timed( Fun(func, B, D; solver=solver) )
                                F = @timed( Fun(func, B, D; solver=solver) )

                                @printf("%3.2e s\t %3.2e bytes",F[2],F[3])
                                @test  msqerror_tol(func, F[1], tol=sqrt(eps(ELT))*10)
                                if func == f
                                    print("\t\t")
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    test_bigfloat = false
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
                    for (func,CT) in ((f,BigFloat),(g,Complex{BigFloat}))
                        B = Basis(71, -T, T, CT)
                        F = @timed( Fun(func, B, D; solver = FE.FE_DirectSolver) )

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

        for D in [Disk(), Disk(1.2,[-0.1,-0.2]), Cube((-1.0,-1.5),(0.5,0.7))]
            show(D); println()

            for solver in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
                show(solver); println()

                for n in ((11,11),)
                    println("\tN = $n")
                    for T in (Extensive ? (FE.default_frame_T(D, Basis),) : ((1.7,1.7),FE.default_frame_T(D, Basis),(2.3,2.3)))
                        print("T = $T\t")
                        B = Basis(n[1],-T[1],T[1]) ⊗ Basis(n[2],-T[2],T[2])
                        for func in (f,g)
                            F = @timed( Fun(func, B, D; solver=solver))

                            @printf("%3.2e s\t %3.2e bytes",F[2],F[3])
                            @test msqerror_tol(func, F[1], tol=1e-3)
                            if func==f
                                print("\t\t")
                                B = Basis(n[1],-T[1],T[1],Complex{Float64}) ⊗ Basis(n[2],-T[2],T[2],Complex{Float64})
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

        for D in (Cube((-0.2,-1.0,0.0),(1.0,-0.1,1.2)),FE.TensorProductDomain(Interval(-1.0,1.0),Disk(1.05)), FE.Ball(1.2,[-0.3,0.25,0.1]))
            show(D); println()
            for solver in (FE.FE_ProjectionSolver, )
                show(solver); println()

                n = FE.default_frame_n(D, Basis)
                println("\tN = $n")

                for T in ((1.7,1.7,1.7), FE.default_frame_T(D, Basis))
                    print("T = $T\t")
                    B = Basis(n[1],-T[1],T[1]) ⊗ Basis(n[2],-T[2],T[2]) ⊗ Basis(n[3],-T[3],T[3])
                    F = @timed( Fun(f, B, D; solver=solver, cutoff=10.0^(3/4*log10(eps(numtype(B))))))
                    @printf("%3.2e s\t %3.2e bytes", F[2], F[3])
                    @test msqerror_tol(f, F[1], tol=1e-2)
                end
            end
        end
    end

    delimit("Random circles")
    dom = FE.randomcircles(10)
    b = FourierBasis(21) ⊗ FourierBasis(21)
    f(x,y) = cos(20*x+22*y)
    @time F = Fun(f,b,dom)
end

# Diagnostics
println()
println("Succes rate:\t$successes/$(successes+failures+errors)")
println("Failure rate:\t$failures/$(successes+failures+errors)")
println("Error rate:\t$errors/$(successes+failures+errors)")

(errors+failures)==0 || error("A total of $(failures+errors) tests failed")

end

