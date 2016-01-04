# A test suite for support functions
module test_suite_support


using BasisFunctions
using FrameFuns
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
global failures=0
global successes=0
global errors=0
# Custom test handler
custom_handler(r::Test.Success) = begin print_with_color(:green, "#\tSuccess "); println("on $(r.expr)"); global successes+=1;  end
custom_handler(r::Test.Failure) = begin print_with_color(:red, "\"\tFailure "); println("on $(r.expr)\""); global failures+=1; end
custom_handler(r::Test.Error) = begin println("\"\t$(typeof(r.err)) in $(r.expr)\""); global errors+=1; end
#custom_handler(r::Test.Error) = Base.showerror(STDOUT,r); 



function delimit(s::AbstractString)
    println("############")
    println("# ",s)
    println("############")
end

#######
# Testing
#######


function test_subgrids()
    delimit("Grid functionality")

    n = 20
    grid1 = EquispacedGrid(n, -1.0, 1.0)
    subgrid1 = FE.MaskedGrid(grid1, Interval(-0.5, 0.7))
    subgrid2 = FE.IndexSubGrid(grid1, 4, 12)

    G1 = EquispacedGrid(n, -1.0, 1.0)
    G2 = EquispacedGrid(n, -1.0, 1.0)
    TensorG = G1 ⊗ G2
    C = Disk(1.0)
    circle_grid = FE.MaskedGrid(TensorG, C)
    @test (length(circle_grid)/length(TensorG)-pi*0.25) < 0.01

    G1s = FE.IndexSubGrid(G1,2,4)
    G2s = FE.IndexSubGrid(G2,3,5)
    TensorGs = G1s ⊗ G2s
    @test G1s[1] == G1[2]
    @test G2s[1] == G2[3]
    @test TensorGs[1,1] == [G1[2],G2[3]]

    # Generic tests for the subgrids
    for (grid,subgrid) in ( (grid1,subgrid1), (grid1,subgrid2), (TensorG, circle_grid))
        print("Subgrid is ")
        println(typeof(subgrid))
        # Count the number of elements in the subgrid
        cnt = 0
        for i in 1:length(grid)
            if i ∈ subgrid
                cnt += 1
            end
        end
        @test cnt == length(subgrid)

        space = DiscreteGridSpace(grid)
        subspace = DiscreteGridSpace(subgrid)
        R = restriction_operator(space, subspace)
        E = extension_operator(subspace, space)

        e = random_expansion(subspace)
        e_ext = E * e
        # Are the elements in the right place?
        cnt = 0
        diff = 0.0
        for i in 1:length(grid)
            if i ∈ subgrid
                cnt += 1
                diff += abs(e[cnt] - e_ext[i])
            end
        end
        @test diff < 1e-6

        e_rest = R * e_ext
        @test sum([abs(e[i]-e_rest[i]) for i in 1:length(e)]) < 1e-6
    end
end




Test.with_handler(custom_handler) do

    test_subgrids()


    delimit("Domains")
    # Integer interval is apparently impossible
    #try Intervala=Interval(-1,1) catch y; message(y) end
    # Interval
    Intervala=Interval(-1.0,1.0)
    Intervala=Intervala+2
    @test FE.left(Intervala)==1
    @test FE.left(2*Intervala)==2
    @test FE.right(Intervala/4)==0.75
    # Disk
    C = Disk(2.0)
    @test FE.in([1.4, 1.4], C)
    @test !FE.in([1.5, 1.5], C)
    @test FE.boundingbox(C) == FE.BBox((-2.0,-2.0),(2.0,2.0))
    # This is certainly unwanted behavior! Due to method inheritance
    @test typeof(1.2*C)==typeof(C*1.2)
    # This is due to a wrong implementation in Scaled Domain
    @test FE.in([1.5,1.5],1.2*C)
    @test FE.in([1.5,1.5],C*1.2)

    #Square
    D = Cube(2)
    @test FE.in([0.9, 0.9],D)
    @test !FE.in([1.1, 1.1],D)
    @test FE.boundingbox(D)==FE.BBox((-1.0,-1.0),(1.0,1.0))
    DS=FE.union(D,C)
    #Cube
    D=Cube((-1.5,0.5,-3.0),(2.2,0.7,-1.0))
    @test FE.in([0.9, 0.6, -2.5],D)
    @test !FE.in([0.0, 0.6, 0.0],D)
    @test FE.boundingbox(D)==FE.BBox((-1.5,0.5,-3.0),(2.2,0.7,-1.0))

    #Ball
    S=Ball(2.0)
    @test FE.in([1.9,0.0,0.0],S)
    @test FE.in([0,-1.9,0.0],S)
    @test FE.in([0.0,0.0,-1.9],S)
    @test !FE.in([1.9,1.9,0.0],S)
    @test FE.boundingbox(S)==FE.BBox((-2.0,-2.0,-2.0),(2.0,2.0,2.0))
    # joint domain
    DS=FE.union(D,S)
    @test FE.in([0.0,0.6,0.0],DS)
    @test FE.in([0.9, 0.6, -2.5],DS)
    @test FE.boundingbox(DS)==FE.BBox((-2.0,-2.0,-3.0),(2.2,2.0,2.0))
    # domain intersection
    DS=FE.intersect(D,S)
    @test !FE.in([0.0,0.6,0.0],DS)
    @test FE.in([0.2, 0.6, -1.1],DS)
    @test FE.boundingbox(DS)==FE.BBox((-1.5,0.5,-2.0),(2.0,0.7,-1.0))
    # domain difference
    DS=D-S
    @test FE.boundingbox(DS)==FE.boundingbox(D)
    # TensorProductDomain 1
    T=FE.tensorproduct(Interval(-1.0,1.0),2)
    FE.in([0.5,0.5],T)
    @test FE.in([0.5,0.5],T)
    @test !FE.in([-1.1,0.3],T)

    # TensorProductDomain 2
    T=FE.TensorProductDomain(Disk(1.05),Interval(-1.0,1.0))
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

    rgrid = FE.IndexSubGrid(grid2, 1, 2*n)
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

    
    
    # Provisional: Restriction transposed is Zeropadding and vice versa
    @test f_restriction'*coef_restricted==f_extension*coef_restricted

    # Provisional:    
    op  = t_restriction * itransform2 * f_extension
    opt = f_restriction * transform2 * t_extension
    @test op*coef_restricted==(opt')*coef_restricted

    ## delimit("Memory Allocation")
    ## delimit("1D")
    ## coef_src = rand(length(grid1))+1im*rand(length(grid1))
    ## coef_dest = rand(length(rgrid))+1im*rand(length(rgrid))
    ## coef2=rand(length(grid2))+1im*rand(length(grid2))

    ## # Until precompilation
    ## FE.apply!(S,coef_src)
    ## FE.apply!(transform2,coef2)
    ## FE.apply!(op,coef_dest,coef_src)
    ## FE.apply!(opt,coef_src,coef_dest)
    ## # In place operators don't allocate (much) memory ()
    ## @test @allocated(FE.apply!(S,coef_src))<100
    ## @test @allocated(FE.apply!(transform2,coef2))<1081
    ## @test @allocated(FE.apply!(op,coef_dest,coef_src)) <1081
    ## @test @allocated(FE.apply!(opt,coef_src,coef_dest)) <1081
    ## delimit("2D")
    ## C=Disk(1.0)
    ## for n=[10,100,200]
    ##     problem = FE.discretize_problem(C,(n,n),(2.0,2.0),(2.0,2.0),FourierBasis,Complex{Float64})
    ##     op=operator(problem)
    ##     opt=FE.operator_transpose(problem)
    ##     coef_src = rand(size(FE.frequency_basis(problem)))+1im*rand(size(FE.frequency_basis(problem)))
    ##     coef_dest = rand(size(FE.time_basis_restricted(problem)))+1im*rand(size(FE.time_basis_restricted(problem)))
    ##     if n==10
    ##         A=FE.matrix(op)
    ##         @test size(A)==size(op)
    ##     end
    ##     FE.apply!(op,coef_dest,coef_src)
    ##     FE.apply!(opt,coef_src,coef_dest)
    ##     mem = @allocated(FE.apply!(op,coef_dest,coef_src))
    ##     println(mem)
    ##     @test @allocated(FE.apply!(op,coef_dest,coef_src)) < 4500
    ##     @test @allocated(FE.apply!(opt,coef_src,coef_dest)) < 4500
    ## end
end

# Diagnostics
println()
println("Succes rate:\t$successes/$(successes+failures+errors)")
println("Failure rate:\t$failures/$(successes+failures+errors)")
println("Error rate:\t$errors/$(successes+failures+errors)")
(errors+failures)==0 || error("A total of $(failures+errors) tests failed")

end

