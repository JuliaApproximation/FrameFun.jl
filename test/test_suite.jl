module test_suite

using DomainSets, BasisFunctions, FrameFun, StaticArrays
using Test, Printf, LinearAlgebra, Random

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
const include_bigfloat_tests = true

########
# Auxiliary functions
########

const v = DomainSets.TypeFactory{SVector}()

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

show_mv_timings(basis, domain, fun; options...) = show_mv_times && show_mv_timings(ExtensionFrame(domain, basis), fun; options...)

function show_mv_timings(dict::Dictionary, fun; options...)
    A, B = FrameFun.discretization(FrameFun.promote_dictionary(dict, fun), fun; options...)
    op = A
    c1 = zeros(eltype(op), size(src(op)))
    c2 = zeros(eltype(op), size(dest(op)))
    apply!(op, c2, c1)
    t = @timed apply!(op, c2, c1)
    if t[3] > 700
        printstyled("\"\tMV product: ", color=:red)
    elseif t[3] == 0
        printstyled("\"\tMV product: ", color=:green)
    else
        printstyled("\"\tMV product: ", color=:blue)
    end
    println(t[3], " bytes allocated, in ", t[2], " seconds")

    global total_mv_allocs += t[3]
    global total_mv_time += t[2]
    return nothing
end



#######
# Testing -- Accuracy
#######

function test_1d_cases()
    delimit("1D")

    ## The tridiagonal solver for Fourier extensions
    solverstyle = TridiagonalProlateStyle()

    n = 61
    T = 1.9
    domain = Interval(-1.0, 1.0)
    basis = FourierBasis(n, -T, T)

    println()
    println("Testing \t solver = $solverstyle, \n\t\t Domain = $domain, \n\t\t Basis = Fourier,\n\t\t ELT = Float64 ")
    verbose && println("N\t T\t Complex?\t abserror\t time\t\t memory   ")

    # Complex and real functions / Float64
    f(x) = cos(x^2-1//2)-1
    g(x) = 1im*cos(x^2-1//2)-1
    for func in (f,g)
        F = @timed Fun(func, basis, domain; solverstyle=solverstyle)
        error = abserror(func, F[1])
        if verbose
            print("$n\t $T\t\t")
            func==g ? print("Y\t") : print("N\t")
            @printf("%3.2e\t %3.2e s\t %3.2e bytes \n",error, F[2],F[3])
        end
        @test  (error < sqrt(eps(Float64))*10)
        show_mv_timings(basis, domain, func; solverstyle=solverstyle)
    end

    @testset "result" for ELT in (Float32,Float64),
            Basis in (FourierBasis, ChebyshevBasis),
            domain in [Interval(ELT(-1.5),ELT(0.7)), DomainSets.UnionDomain(Interval(ELT(-1.5),ELT(-0.5)),Interval(ELT(0.5),ELT(1.5)))],
            solverstyle in (AZStyle(), DirectStyle())

        println()
        println("Testing \t solver = $solverstyle, \n\t\t Domain = $domain, \n\t\t Basis = $(name(instantiate(Basis,10))),\n\t\t ELT = $ELT ")
        verbose && println("N\t T\t Complex?\t abserror\t time\t\t memory   ")

        for n in (61,)

            # There is some symmetry around T=2, test smaller and larger values
            for T in (ELT(1.9),)
                for func in (f,g)

                    basis = Basis{ELT}(n, -T, T)
                    F = @timed Fun(func, basis, domain; solverstyle=solverstyle)
                    error = abserror(func, F[1])
                    if verbose
                        print("$n\t $T\t\t")
                        func==g ? print("Y\t") : print("N\t")
                        @printf("%3.2e\t %3.2e s\t %3.2e bytes \n",error, F[2],F[3])
                    end
                    @test  (error < sqrt(eps(ELT))*10)
                    show_mv_timings(basis, domain, func; solverstyle=solverstyle)
                end
            end
        end
    end
end

function test_bigfloat()
    f(x) = cos(x.^2) - big(1.0)
    g(x) = big(1.0)im * cos(x.^2) - big(1.0)
    @testset "result" for Basis in (FourierBasis, ChebyshevBasis),
            domain in [Interval(big(-3.0)/2, big(7.0)/10)]
        println()
        println("Testing \t solver = QR_solver\n\t\t Domain = $domain, \n\t\t Basis = $(name(instantiate(Basis,10))),\n\t\t ELT = BigFloat ")
        verbose && println("N\t T\t Complex?\t abserror\t time\t\t memory   ")
        for T in (big(17.0)/10,)
            for func in (f,g)
                basis = Basis(91, -T, T)
                F = @timed Fun(func, basis, domain; solverstyle = DirectStyle(), directsolver=:qr)
                error = abserror(func, F[1])
                if verbose
                    @printf("91\t %3.2e\t\t",T)
                    func==g ? print("Y\t") : print("N\t")
                    @printf("%3.2e\t %3.2e s\t %3.2e bytes \n",error, F[2],F[3])
                end
                @test error < 1e-20
                show_mv_timings(basis, domain, func; solverstyle = DirectStyle(), directsolver=:qr)
            end
        end
    end
end

function test_2d_cases()
    delimit("2D")

    f(x,y) = cos(0.5*x)+2*sin(0.2*y)-1.0*x*y
    g(x,y) = 1im*cos(0.5*x)+2*sin(0.2*y)-1.0im*x*y
    @testset "result" for Basis in (FourierBasis, ChebyshevBasis),
            domain in [disk(1.2,v[-0.1,-0.2]), cube((-1.0,-1.5),(0.5,0.7))],
            solverstyle in (AZStyle(), DirectStyle())

        println()
        println("Testing \t solver = $solverstyle \n\t\t Domain = $domain, \n\t\t Basis = $(name(instantiate(Basis,10)⊗instantiate(Basis,10))),\n\t\t ELT = Float64 ")
        verbose && println("N\t\t T\t\t Complex?\t abserror\t time\t\t memory   ")

        for n in ((11,11),)
            for T in ((1.7,2.3),)
                basis = Basis(n[1],-T[1],T[1]) ⊗ Basis(n[2],-T[2],T[2])
                for func in (f,g)
                    F = @timed Fun(func, basis, domain; solverstyle=solverstyle)
                    error = abserror(func, F[1])
                    if verbose
                        print("$n \t $T \t\t")
                        func==g ? print("Y\t") : print("N\t")
                        @printf("%3.2e\t %3.2e s\t %3.2e bytes \n",error, F[2],F[3])
                    end
                    @test error < 1e-3
                    show_mv_timings(basis, domain, func; solverstyle=solverstyle)
                end

            end
        end
    end
end

function test_3d_cases()
    delimit("3D")

    f(x,y,z) = cos(x)+sin(y)-x*z
    # @testset "result" for Basis in (FourierBasis, ChebyshevBasis), D in (Cube((-1.2,-1.0,-0.9),(1.0,0.9,1.2)),FrameFun.tensorproduct(Interval(-1.0,1.0),Disk(1.05)), FrameFun.Ball(1.2,[-0.3,0.25,0.1])), solver in (AZSolver, )
    #             show(solver); println()
    @testset "result" for Basis in (FourierBasis, ChebyshevBasis),
                domain in (cube((-1.2,-1.0,-0.9),(1.0,0.9,1.2)), ball(1.2,v[-0.3,0.25,0.1])),
                solverstyle in (AZStyle(), )
        show(solverstyle); println()
        println()
        println("Testing \t solver = $solverstyle \n\t\t Domain = $domain, \n\t\t Basis = $(name(instantiate(Basis,10)⊗instantiate(Basis,10)⊗instantiate(Basis,10))),\n\t\t ELT = Float64 ")
        verbose && println("N\t\t T\t\t Complex?\t abserror\t time\t\t memory   ")

        n = (12,12,12)

        for T in ((1.7,1.2,1.3),)
            basis = Basis(n[1],-T[1],T[1]) ⊗ Basis(n[2],-T[2],T[2]) ⊗ Basis(n[3],-T[3],T[3])
            F = @timed Fun(f, basis, domain; solverstyle=solverstyle, threshold=10.0^(3/4*log10(eps(real(codomaintype(basis))))), oversamplingfactor=1.5)
            error = abserror(f, F[1])
            if verbose
                print("$n \t $T \t\t")
                print("N\t")
                @printf("%3.2e\t %3.2e s\t %3.2e bytes \n",error, F[2],F[3])
            end
            @test  error < 1e-1
            show_mv_timings(basis, domain, f; solverstyle=solverstyle, threshold=10.0^(3/4*log10(eps(real(codomaintype(basis))))), oversamplingfactor=1.5)
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

# delimit("Random circles")
# dom = FrameFun.randomcircles(10)
# #    b = FourierBasis(21) ⊗ FourierBasis(21)
# b = FourierBasis(31) ⊗ ChebyshevBasis(31)
# f(x,y) = cos(20*x+22*y)
# @time F = Fun(f,b,dom,verbose=true)

if show_mv_times
    println("Total bytes in MV products:\t$total_mv_allocs")
    println("Total time in MV products:\t$total_mv_time")
end


end
