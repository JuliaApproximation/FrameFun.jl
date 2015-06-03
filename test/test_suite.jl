using BasisFunctions
using FrameFuns
using FactCheck
# Parameters
show_backtrace = false

# Algorithm accuracy below tol
function msqerror{N,T}(f::Function,F::Fun{N,T},vals=100,tol=1e-12)
    for i in range(1,vals)
        x=rand()

facts("Algorithm implementation and accuracy") do

    context("1D") do

        f(x)=x
        g(x)=x+1im*x

        # Standard methods
        F=ExpFun(f)
        @fact msqerror_tol(f,F)
        G=ExpFun(g,n=100,T=1.1)
        @fact msqerror_tol(g,G)

    end

    context("2D") do

        f(x)=x[1]+2*x[2]
        g(x)=x[1]-2im*x[2]

        # Standard methods -- Default domain results in BoundsErrors
        try
            F = ExpFun(f)
            G = ExpFun(g,n=20,)
            @fact msqerror_tol(f,F)=>less_than(1e-8)
            @fact msqerror_tol(f,F)=>less_than(1e-8)

        catch y
            println(y)
            show_backtrace && Base.show_backtrace(STDOUT,catch_backtrace())
        end
            

        # Standard methods -- Circle
        C = Circle(1.0)
        F=ExpFun(f,C)
        @fact msqerror_tol(f,F)=>less_than(1e-8)
        G=ExpFun(g,n=20,T=1.4,C)
        @fact msqerror_tol(g,G)=>less_than(1e-8)
    end
end

facts("Operator functionality") do
    T=Float64
    n=21
    m=41
    l=80
    N=1
    t = (l[1]*one(T)) / ((m[1]-1)*one(T))

    fbasis1 = FourierBasis(n[1], -one(T), one(T) + 2*(t-1))
    fbasis2 = FourierBasis(l[1], -one(T), one(T) + 2*(t-1))
    tens_fbasis1 = tensorproduct(fbasis1, N)
    tens_fbasis2 = tensorproduct(fbasis2, N)

    grid1 = grid(fbasis1)
    grid2 = grid(fbasis2)
    tens_grid1 = tensorproduct(grid1, N)
    tens_grid2 = tensorproduct(grid2, N)

    rgrid = FrameFuns.EquispacedSubGrid(grid2, 1, m[1])
    tens_rgrid = tensorproduct(rgrid, N)

    tbasis1 = TimeDomain(grid1)
    tbasis2 = TimeDomain(grid2)
    tens_tbasis1 = tensorproduct(tbasis1, N)
    tens_tbasis2 = tensorproduct(tbasis2, N)

    tbasis_restricted = TimeDomain(rgrid)
    tens_tbasis_restricted = tensorproduct(tbasis_restricted, N)

    f_extension = ZeroPadding(tens_fbasis1, tens_fbasis2)
    f_restriction = Restriction(tens_fbasis2, tens_fbasis1)

    t_extension = ZeroPadding(tens_tbasis_restricted, tens_tbasis2)
    t_restriction = Restriction(tens_tbasis2, tens_tbasis_restricted)

    transform1 = FastFourierTransformFFTW(tens_tbasis1, tens_fbasis1)
    itransform1 = InverseFastFourierTransformFFTW(tens_fbasis1, tens_tbasis1)

    transform2 = FastFourierTransformFFTW(tens_tbasis2, tens_fbasis2)
    itransform2 = InverseFastFourierTransformFFTW(tens_fbasis2, tens_tbasis2)

end




