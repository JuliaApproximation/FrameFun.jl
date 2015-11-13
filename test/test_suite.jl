module test_suite


using BasisFunctions
using FrameFuns
using Base.Test
FE=FrameFuns
BA=BasisFunctions
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
function msqerror_tol{N,T}(f::Function,F::FE.Fun{N,T};vals::Int=200,tol=1e-6)
    # Find the closest bounding grid around the domain
    TB=FE.box(FE.domain(F))
    
    point=Array{T}(N)
    elements=0
    error=0
    l=left(TB)
    r=right(TB)
    
    for i in 1:vals
        point=l+(r-l).*rand(T,N)
        if FE.in(point,FE.domain(F))
            elements+=1
            error+=abs(f(point)-F(point...))
        end
    end
    @printf(" %3.2e",error/elements)
    return error>0 ? error/elements<tol : false
end
# I don't like these functions, replace by boundingbox?
boundinggrid(g::AbstractGrid) = g
boundinggrid(g::FE.MaskedGrid) = g.grid

function message(y)
    println("\"\t",typeof(y),"\"")
    Base.showerror(STDOUT,y)
    global errors+=1
end

function message(y, backtrace)
    rethrow(y)
    print("\tFailure due to ",typeof(y))
    Base.showerror(STDOUT,y)
    Base.show_backtrace(STDOUT,backtrace)
    global errors+=1
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
    delimit("grid functionality")

    delimit("MaskedGrid")
    G1=BA.EquispacedGrid(100,-1.0,1.0)
    G2=BA.EquispacedGrid(100,-1.0,1.0)
    TensorG=BA.TensorProductGrid(G1,G2)
    C=Circle(1.0)
    G4=FE.MaskedGrid(TensorG,C)
    @test (length(G4)/length(TensorG)-pi*0.25)<0.01
    # I'm assuming here MaskedGrids aren't supposed to be indexed.
    @test_throws Exception G4[1,1]

    delimit("SubGrid")
    G1s=FE.IndexedSubGrid(G1,2,4)
    G2s=FE.IndexedSubGrid(G2,3,5)
    @test G1s[1]==G1[2]
    @test G2s[1]==G2[3]
    TensorGs=TensorProductGrid(G1s,G2s)
    @test TensorGs[1,1]==[G1[2],G2[3]]

    delimit("Domains")
    # Integer interval is apparently impossible
    #try Intervala=Interval(-1,1) catch y; message(y) end
    # Interval
    Intervala=Interval(-1.0,1.0)
    Intervala=Intervala+2
    @test FE.left(Intervala)==1
    @test FE.left(2*Intervala)==2
    @test FE.right(Intervala/4)==0.75
    # Circle
    C=Circle(2.0)
    @test FE.in([1.4, 1.4],C)
    @test !FE.in([1.5, 1.5],C)
    @test FE.box(C)==FE.BBox((-2.0,-2.0),(2.0,2.0))
    # This is certainly unwanted behavior! Due to method inheritance
    @test typeof(1.2*C)==typeof(C*1.2)
    # This is due to a wrong implementation in Scaled Domain
    @test FE.in([1.5,1.5],1.2*C)
    @test FE.in([1.5,1.5],C*1.2)
    #Square
    D=Cube(2)
    @test FE.in([0.9, 0.9],D)
    @test !FE.in([1.1, 1.1],D)
    @test FE.box(D)==FE.BBox((-1.0,-1.0),(1.0,1.0))
    DS=FE.join(D,C)
    #Cube
    D=Cube((-1.5,0.5,-3.0),(2.2,0.7,-1.0))
    @test FE.in([0.9, 0.6, -2.5],D)
    @test !FE.in([0.0, 0.6, 0.0],D)
    @test FE.box(D)==FE.BBox((-1.5,0.5,-3.0),(2.2,0.7,-1.0))

    #Sphere
    S=Sphere(2.0)
    @test FE.in([1.9,0.0,0.0],S)
    @test FE.in([0,-1.9,0.0],S)
    @test FE.in([0.0,0.0,-1.9],S)
    @test !FE.in([1.9,1.9,0.0],S)
    @test FE.box(S)==FE.BBox((-2.0,-2.0,-2.0),(2.0,2.0,2.0))
    # joint domain
    DS=FE.join(D,S)
    @test FE.in([0.0,0.6,0.0],DS)
    @test FE.in([0.9, 0.6, -2.5],DS)
    @test FE.box(DS)==FE.BBox((-2.0,-2.0,-3.0),(2.2,2.0,2.0))
    # domain intersection
    DS=FE.intersect(D,S)
    @test !FE.in([0.0,0.6,0.0],DS)
    @test FE.in([0.2, 0.6, -1.1],DS)
    @test FE.box(DS)==FE.BBox((-1.5,0.5,-2.0),(2.0,0.7,-1.0))
    # domain difference
    DS=D-S
    @test FE.box(DS)==FE.box(D)
    # TensorProductDomain 1
    T=FE.tensorproduct(Interval(-1.0,1.0),2)
    FE.in([0.5,0.5],T)
    @test FE.in([0.5,0.5],T)
    @test !FE.in([-1.1,0.3],T)

    # TensorProductDomain 2
    T=FE.TensorProductDomain(Circle(1.05),Interval(-1.0,1.0))
    @test FE.in([0.5,0.5,0.8],T)
    @test !FE.in([-1.1,0.3,0.1],T)

    delimit("Basis and operator functionality")

    
    n=100
    a=-1.2
    b=0.7
    delimit("Fourier Basis")
    fbasis1 = FourierBasis(n+1,a,b)
    @test grid(fbasis1)==BA.PeriodicEquispacedGrid(n+1,a,b)
    r=rand()
    @test abs(fbasis1(2,r)-exp(2*pi*1im*(2-1)*(r-a)/(b-a)))<1e-13
    r=sqrt(BigFloat(pi))
    @test abs(fbasis1(2,r)-exp(2*pi*1im*(2-1)*(r-a)/(b-a)))<1e-13
    a=BigFloat(-1.0); b=BigFloat(1.0)
    fbasis1 = FourierBasis(n+1,a,b)
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

    rgrid = FE.IndexedSubGrid(grid2, 1, 2*n)
    tbasis1 = DiscreteGridSpace(grid1)
    tbasis2 = DiscreteGridSpace(grid2)

    tbasis_restricted = DiscreteGridSpace(rgrid)

    f_extension = Extension(fbasis1, fbasis2)
    f_restriction = Restriction(fbasis2, fbasis1)

    t_extension = Extension(tbasis_restricted, tbasis2)
    t_restriction = Restriction(tbasis2, tbasis_restricted)
   
    transform1 = transform_operator(tbasis1, fbasis1)
    itransform1 = transform_operator(fbasis1, tbasis1)

    transform2 = transform_operator(tbasis2, fbasis2)
    itransform2 = transform_operator(fbasis2, tbasis2)

    I=IdentityOperator(tbasis1)
    S=ScalingOperator(tbasis1,2)

    coef_src=rand(length(grid1))
    coef2=rand(length(grid2))

    @test I*coef_src==coef_src
    @test S*coef_src==2*coef_src
    # Verify transposes of Identity and Scaling operators
    @test I'*coef_src==coef_src
    @test S'*coef_src==2*coef_src 

    # Extension and Restriction nullify (this doesn't work for even n!)
    coef_restricted = rand(length(rgrid))
    @test t_restriction*t_extension*coef_restricted==coef_restricted
    coef_restricted = rand(length(grid1))+1im*rand(length(grid1))
    @test f_restriction*f_extension*coef_restricted==coef_restricted
    @test (f_restriction*f_extension)*coef_restricted==f_restriction*(f_extension*coef_restricted)

    # Trying to do a Fourier Transform on real numbers doesn't work because of Typing.
    @test FE.apply!(transform2,coef2) 
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
    FE.apply!(S,coef_src)
    FE.apply!(transform2,coef2)
    FE.apply!(op,coef_dest,coef_src)
    FE.apply!(opt,coef_src,coef_dest)
    # In place operators don't allocate (much) memory ()
    @test @allocated(FE.apply!(S,coef_src))<100
    @test @allocated(FE.apply!(transform2,coef2))<881
    @test @allocated(FE.apply!(op,coef_dest,coef_src)) <881
    @test @allocated(FE.apply!(opt,coef_src,coef_dest)) <881
    delimit("2D")
    C=Circle(1.0)
    for n=[10,100,200]
        problem = FE.discretize_problem(C,(n,n),(2.0,2.0),(2.0,2.0),FourierBasis,Complex{Float64})
        op=operator(problem)
        opt=FE.operator_transpose(problem)
        coef_src = rand(size(FE.frequency_basis(problem)))+1im*rand(size(FE.frequency_basis(problem)))
        coef_dest = rand(size(FE.time_basis_restricted(problem)))+1im*rand(size(FE.time_basis_restricted(problem)))
        if n==10
            A=FE.matrix(op)
            @test size(A)==size(op)
        end
        FE.apply!(op,coef_dest,coef_src)
        FE.apply!(opt,coef_src,coef_dest)
        mem = @allocated(FE.apply!(op,coef_dest,coef_src))
        println(mem)
        @test @allocated(FE.apply!(op,coef_dest,coef_src)) < 4500
        @test @allocated(FE.apply!(opt,coef_src,coef_dest)) < 4500
    end

    
    delimit("Algorithm Implementation and Accuracy")
    delimit("1D")
    
    f(x)=x[1]-1.0
   for funtype in (ExpFun, ChebyFun)
        println("Fun Type: ",funtype)
        for D in [FE.default_fourier_domain_1d() Interval(-1.5,0.7) Interval(-1.5,-0.5)+Interval(0.5,1.5)]
            show(D); print("\n")
            for solver_type in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
                show(solver_type);print("\n")
                for n in [FE.default_fourier_n(D) 49]
                    println("\tN = $n")
                    for T in [1.7 FE.default_fourier_T(D) 2.3]
                        print("T = $T \t")
                        try
                            F=@timed(funtype(f,D,solver_type,n=n,T=T))
                            @printf("%3.2e s\t %3.2e bytes",F[2],F[3])
                            @test  msqerror_tol(f,F[1],tol=1e-7)
                        catch y
                            message(y)
                        end
                    end
                end
            end
        end
    end
    delimit("2D") 

    f(x)=x[1]+2*x[2]-1.0
    f(x,y)=x+2*y-1.0
#    for funtype in (ExpFun,ChebyFun)
    for funtype in (ExpFun,)
        println("Fun Type: ",funtype)
        for D in [Cube((-1.0,-1.5),(0.5,0.7)) Circle(1.0) Circle(2.0,[-2.0,-2.0])]       
            show(D); print("\n")
            for solver_type in (FE.FE_DirectSolver, FE.FE_ProjectionSolver)
                show(solver_type);print("\n")
                for n in (FE.default_fourier_n(D),(12,12))
                    println("\tN = $n")
                    for T in ((1.7,1.7),FE.default_fourier_T(D),(2.3,2.3))
                        print("T = $T\t")
                        try
                            F=@timed(funtype(f,D,solver_type,n=n,T=T))
                            @printf("%3.2e s\t %3.2e bytes",F[2],F[3])
                            @test msqerror_tol(f,F[1],tol=1e-5)
                        catch y
                            message(y,catch_backtrace())
                        end
                    end
                end
            end
        end
    end
    return
    delimit("3D")

    f(x)=x[1]+x[2]-x[3]
    
    for D in (FE.TensorProductDomain(Interval(-1.0,1.0),Circle(1.05)),Cube((-0.2,-0.3,0.0),(1.0,-0.1,1.2)),Cylinder(1.05,1.2), FE.Sphere(1.2,[-1.3,0.25,1.0]))        
        show(D); print("\n")
        for solver_type in (FE.FE_ProjectionSolver, FE.FE_DirectSolver)
            show(solver_type);print("\n")
            for n in ((3,3,3), FE.default_fourier_n(D))
                println("\tN = $n")
                for T in ((1.7,1.7,1.7), FE.default_fourier_T(D), (2.3,2.3,2.3))
                    print("T = $T\t")
                    try
                        F=@timed(ExpFun(f,D,solver_type,n=n,T=T))
                        @printf("%3.2e s\t %3.2e bytes",F[2],F[3])
                        @test msqerror_tol(f,F[1],tol=1e-2)
                    catch y
                        message(y)
                    end
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
end

