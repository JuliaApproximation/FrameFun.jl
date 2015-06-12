module test_suite

using BasisFunctions
using FrameFuns
using Base.Test

########
# Auxiliary functions
########

# Keep track of successes, failures and errors
global failures=0
global successes=0
global errors=0
# Custom test handler
custom_handler(r::Test.Success) = begin println("#\tSucces on $(r.expr)"); global successes+=1;  end
custom_handler(r::Test.Failure) = begin println("\"\tFailure on $(r.expr)\""); global failures+=1; end
custom_handler(r::Test.Error) = begin println("\"\t$(typeof(r.err)) in $(r.expr)\""); global errors+=1; end
#custom_handler(r::Test.Error) = Base.showerror(STDOUT,r); 


# Algorithm accuracy below tol
function msqerror_tol{N,T}(f::Function,F::FrameFuns.Fun{N,T};vals::Int=200,tol=1e-6)
    # Find the closest bounding grid around the domain
    g=boundinggrid(grid(FrameFuns.time_basis_restricted(FrameFuns.problem(F))))
    s=size(g)
    x=Array{Int}(N)
    point=Array{T}(N)
    elements=0
    error=0
    # Generate some points inside the domain, and compare with the target function
    for i in 1:vals
        for j in 1:N
            x[j]=rand(1:s[j])
        end
        point = FrameFuns.getindex(g,x...)
        if FrameFuns.in(point,FrameFuns.domain(F))
            elements+=1
            error+=abs(f(point)-F(point...))
        end
    end
    if !(error/elements < tol) @printf(" %3.2e/%3.2e",error/elements,tol) end
    return error>0 ? error/elements<tol : false
end
# I don't like these functions, replace by boundingbox?
boundinggrid(g::AbstractGrid) = g
boundinggrid(g::FrameFuns.MaskedGrid) = g.grid

function message(y)
    println("\"\t",typeof(y),"\"")
    Base.showerror(STDOUT,y)
    global errors+=1
end

function message(y, backtrace)
    rethrow(y)
    println("\tFailure due to ",typeof(y))
    Base.showerror(STDOUT,y)
    Base.show_backtrace(STDOUT,backtrace)
    global errors+=1
end

function delimit(s::String)
    println("############")
    println("# ",s)
    println("############")
end

#######
# Testing
#######


Test.with_handler(custom_handler) do
    delimit("grid functionality")

    delimit("1D")

    # Careful with these values, grid evaluation is not optimal arithmetically (for performance reasons)
    a=-1.2
    b=0.7
    n=100;
    G1=BasisFunctions.EquispacedGrid(n,a,b)
    G2=BasisFunctions.PeriodicEquispacedGrid(n,a,b)
    for G in [G1 G2]
        @test G[1]==a
        iseven(n) && @test G[round(Int,n/2+1)]==(a+b)/2 
    end
    @test G1[n+1]==b
    @test_throws BoundsError G2[n+1]==b
    @test length(G1)==n+1
    @test length(G2)==n

    delimit("2D")
    # Two different grid types cannot be combined, same type is possible (is this necessary?)
    a2=-1.0
    b2=3.0
    n2=50;
    G3=BasisFunctions.EquispacedGrid(n2,a2,b2)
    for G in [G2 G3]
        try
            TensorG=TensorProductGrid((G1,G))
            @test length(TensorG)==length(G1)*length(G)
            @test TensorG[1,1]==[a, a2]
            (iseven(n) && iseven(n2)) && @test TensorG[round(Int,n/2+1),round(Int,n2/2+1)]==[(a+b)/2,(a2+b2)/2]
        catch y;
            message(y)
        end
    end

    delimit("MaskedGrid")
    G1=BasisFunctions.EquispacedGrid(100,-1.0,1.0)
    G2=BasisFunctions.EquispacedGrid(100,-1.0,1.0)
    TensorG=BasisFunctions.TensorProductGrid((G1,G2))
    C=Circle(1.0)
    G4=FrameFuns.MaskedGrid(TensorG,C)
    @test (length(G4)/length(TensorG)-pi*0.25)<0.01
    # I'm assuming here MaskedGrids aren't supposed to be indexed.
    @test_throws Exception G4[1,1]

    delimit("SubGrid")
    G1s=FrameFuns.EquispacedSubGrid(G1,2,4)
    G2s=FrameFuns.EquispacedSubGrid(G2,3,5)
    @test G1s[1]==G1[2]
    @test G2s[1]==G2[3]
    TensorGs=TensorProductGrid((G1s,G2s))
    @test TensorGs[1,1]==[G1[2],G2[3]]

    delimit("Domains")
    # Integer interval is apparently impossible
    #try Intervala=Interval(-1,1) catch y; message(y) end
    # Interval
    Intervala=Interval(-1.0,1.0)
    Intervala=Intervala+2
    @test FrameFuns.left(Intervala)==1
    @test FrameFuns.left(2*Intervala)==2
    @test FrameFuns.right(Intervala/4)==0.75
    # Circle
    C=Circle(2.0)
    @test FrameFuns.in([1.4, 1.4],C)
    @test !FrameFuns.in([1.5, 1.5],C)
    # This is certainly unwanted behavior! Due to method inheritance
    @test typeof(1.2*C)==typeof(C*1.2)
    # This is due to a wrong implementation in Scaled Domain
    @test FrameFuns.in([1.5,1.5],1.2*C)
    @test FrameFuns.in([1.5,1.5],C*1.2)
    
    delimit("Basis and operator functionality")

    delimit("Fourier Basis")
    fbasis1 = FourierBasis(n,a,b)
    @test grid(fbasis1)==BasisFunctions.PeriodicEquispacedGrid(n,a,b)
    r=rand()
    @test abs(fbasis1(2,r)-exp(2*pi*1im*(2-1)*(r-a)/(b-a)))<1e-13
    r=sqrt(BigFloat(pi))
    @test abs(fbasis1(2,r)-exp(2*pi*1im*(2-1)*(r-a)/(b-a)))<1e-13
    a=BigFloat(-1.0); b=BigFloat(1.0)
    fbasis1 = FourierBasis(n,a,b)
    r=rand()
    @test abs(fbasis1(2,r)-exp(2*pi*1im*(2-1)*(r-a)/(b-a)))<1e-13
    r=sqrt(BigFloat(pi))
    @test abs(fbasis1(2,r)-exp(2*BigFloat(pi)*1im*(2-1)*(r-a)/(b-a)))<1e-30
    delimit("Operator Functionality")

    
    # This is probably too much operators! Close to half of these should be replaceable by transposes

    a=-1.0
    b=1.0
    n=1001
    fbasis1 = FourierBasis(n, a, b)
    fbasis2 = FourierBasis(4*n, a, b)

    grid1 = grid(fbasis1)
    grid2 = grid(fbasis2)

    rgrid = FrameFuns.EquispacedSubGrid(grid2, 1, 2*n)
    tbasis1 = TimeDomain(grid1)
    tbasis2 = TimeDomain(grid2)

    tbasis_restricted = TimeDomain(rgrid)

    f_extension = ZeroPadding(fbasis1, fbasis2)
    f_restriction = Restriction(fbasis2, fbasis1)

    t_extension = ZeroPadding(tbasis_restricted, tbasis2)
    t_restriction = Restriction(tbasis2, tbasis_restricted)
   
    transform1 = transform_operator(tbasis1, fbasis1)
    itransform1 = transform_operator(fbasis1, tbasis1)

    transform2 = transform_operator(tbasis2, fbasis2)
    itransform2 = transform_operator(fbasis2, tbasis2)

    I=IdentityOperator(tbasis1)
    S=ScalingOperator(2,tbasis1,tbasis1)

    coef_src=rand(length(grid1))
    coef2=rand(length(grid2))

    @test I*coef_src==coef_src
    @test S*coef_src==2*coef_src
    # Identity and Scaling operators don't have transposes
    @test I'*coef_src==coef_src
    @test S'*coef_src==2*coef_src 

    # ZeroPadding and Restriction nullify (this doesn't work for even n!)
    coef_restricted = rand(length(rgrid))
    @test t_restriction*t_extension*coef_restricted==coef_restricted
    coef_restricted = rand(length(grid1))+1im*rand(length(grid1))
    @test f_restriction*f_extension*coef_restricted==coef_restricted
    @test (f_restriction*f_extension)*coef_restricted==f_restriction*(f_extension*coef_restricted)

    # Trying to do a Fourier Transform on real numbers doesn't work because of Typing.
    @test FrameFuns.apply!(transform2,coef2) 
    coef2=rand(length(grid2))+1im*rand(length(grid2))
    # Provisional: Restriction transposed is Zeropadding and vice versa
    @test f_restriction'*coef_restricted==f_extension*coef_restricted

    # Provisional:    
    op  = t_restriction * itransform2 * f_extension
    opt = f_restriction * transform2 * t_extension
    @test op*coef_restricted==opt'*coef_restricted

    delimit("Memory Allocation")
    delimit("1D")
    coef_src = rand(length(grid1))+1im*rand(length(grid1))
    coef_dest = rand(length(rgrid))+1im*rand(length(rgrid))
    # Until precompilation
    FrameFuns.apply!(S,coef_src)
    FrameFuns.apply!(transform2,coef2)
    FrameFuns.apply!(op,coef_dest,coef_src)
    FrameFuns.apply!(opt,coef_src,coef_dest)
    # In place operators don't allocate (much) memory ()
    @test @allocated(FrameFuns.apply!(S,coef_src))<100
    @test @allocated(FrameFuns.apply!(transform2,coef2))<881
    @test @allocated(FrameFuns.apply!(op,coef_dest,coef_src)) <881
    @test @allocated(FrameFuns.apply!(opt,coef_src,coef_dest)) <881
    delimit("2D")
    C=Circle(1.0)
    for n=[10,100,200]
        problem = FrameFuns.default_fourier_problem(C,n,2.0,2.0)
        op=operator(problem)
        opt=FrameFuns.operator_transpose(problem)
        coef_src = rand(size(FrameFuns.frequency_basis(problem)))+1im*rand(size(FrameFuns.frequency_basis(problem)))
        coef_dest = rand(size(FrameFuns.time_basis_restricted(problem)))+1im*rand(size(FrameFuns.time_basis_restricted(problem)))
        println(size(coef_src),size(op),size(coef_dest))
        if n==10
            A=FrameFuns.matrix(op)
            @test size(A)==size(op)
        end
        FrameFuns.apply!(op,coef_dest,coef_src)
        FrameFuns.apply!(opt,coef_src,coef_dest)
        @test @allocated(FrameFuns.apply!(op,coef_dest,coef_src)) <3745
        @test @allocated(FrameFuns.apply!(opt,coef_src,coef_dest)) <3745
    end
    delimit("Algorithm Implementation and Accuracy")
    delimit("1D")

    f(x)=x
    g(x)=x+1im*x

    # The fact that n is set to 2n+1 in default_fourier_problem is a bit worrisome to me.
    for D in [FrameFuns.default_fourier_domain_1d() Interval(-1.5,0.7)]        
        show(D); print("\n")
        for solver_type in (FrameFuns.FE_ProjectionSolver, FrameFuns.FE_DirectSolver, FrameFuns.FE_IterativeSolverLSQR)
            show(solver_type);print("\n")
            for n in [FrameFuns.default_fourier_n(D) 49]
                println("\tN = $n")
                for T in [1.4 FrameFuns.default_fourier_T(D) 2.3]
                    print("T = $T \t")
                    try
                        F=ExpFun(f,D,solver_type,n=n,T=T)
                        @test msqerror_tol(f,F,tol=1e-7)
                    catch y
                        message(y)
                    end
                end
            end
        end
    end
    delimit("2D") 

    f(x)=x[1]+2*x[2]
    
    # Standard methods -- Default domain results in BoundsErrors because it's 1d - Talked about this with Daan
    try
        F = ExpFun(f)
        @test msqerror_tol(f,F)
    catch y
        message(y)
    end
    for D in [Circle(1.0) Circle(3.0) 1.0*Cube(2)]        
        show(D); print("\n")
        for solver_type in (FrameFuns.FE_DirectSolver, FrameFuns.FE_ProjectionSolver, FrameFuns.FE_IterativeSolverLSQR)
            show(solver_type);print("\n")
            for n in [FrameFuns.default_fourier_n(D) 12]
                println("\tN = $n")
                for T in [1.4 FrameFuns.default_fourier_T(D) 2.3]
                    print("T = $T\t")
                    try
                        F=ExpFun(f,D,solver_type,n=n,T=T)
                        @test msqerror_tol(f,F,tol=1e-5)
                    catch y
                        message(y,catch_backtrace())
                    end
                end
            end
        end
    end
    return
    # Circle Test fails because some points that "should" be in the circle have no corresponding point in mask, since the domain is not adapted
    F=ExpFun(f,Circle(2.0))
    @test abs(F(-1.4,0.0)-f([-1.4,0.0]))<1e-1
    # Cube Test fails because of unknown reasons

    delimit("3D")

    f(x)=x[1]+x[2]-x[3]
    for D in [FrameFuns.Sphere() Sphere(3.0) Cube(3)]        
        show(D); print("\n")
        for solver_type in (FrameFuns.FE_ProjectionSolver, FrameFuns.FE_DirectSolver, FrameFuns.FE_IterativeSolverLSQR)
            show(solver_type);print("\n")
            for n in [3 4]
                println("\tN = $n")
                for T in [1.4 FrameFuns.default_fourier_T(D) 2.3]
                    print("T = $T\t")
                    try
                        F=ExpFun(f,D,solver_type,n=n,T=T)
                        @test msqerror_tol(f,F,tol=1e-2)
                    catch y
                        message(y)
                    end
                end
            end
        end
    end

end

# Diagnostics
println("Succes rate:\t$successes/$(successes+failures+errors)")
println("Failure rate:\t$failures/$(successes+failures+errors)")
println("Error rate:\t$errors/$(successes+failures+errors)")
end

