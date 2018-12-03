module test_suite_applications

using DomainSets, BasisFunctions, FrameFun, Plots
# Plots loads the default backend (PyPlot unless otherwise specified)

using Test

FE = FrameFun
BA = BasisFunctions

function delimit(s::AbstractString)
    println()
    println("############")
    println("# ",s)
    println("############")
end

delimit("Plotting")

function test_plots()
    gr()
    # test plotting functionsets
    B = ChebyshevBasis(5)
    plot(B)

    plot(B[4])

    F = FourierBasis(5)
    plot(F[2],plot_complex=true)

    B = ChebyshevBasis(5,-2,0)⊕FourierBasis(5,0,2)
    plot(B[2:8],plot_complex=true)

    F = FourierBasis(10)⊗FourierBasis(10)
    plot(F[82],plot_complex=false)

    G = BasisFunctions.grid(ChebyshevBasis(51))
    plot(G)

    G = BasisFunctions.grid(ChebyshevBasis(51,-1,0.3)⊗ChebyshevBasis(51,-0.5,0.5))
    G = FrameFun.subgrid(G,mandelbrot())
    plot(G)

    B = FourierBasis(21,-1,1)⊗FourierBasis(21,-1,1)⊗FourierBasis(21,-1,1)
    D = ball()\FrameFun.cube(-0.5,0.5,-0.5,0.5,-2.,2.)
    Df = ExtensionFrame(D,B)
    G = BasisFunctions.grid(Df)
    plot(G,size=(400,400))

    basis = FourierBasis(51,-2.0,-0.5)
    domain = -1.7..1.0
    f(x) = cos(3*x.^2)
    F = Fun(f, basis, domain)
    # Easily combine multiple plots
    plot(BasisFunctions.grid(dictionary(F)),label="grid",markercolor=:white)
    plot!(F,label="function", plot_ext=true)
    plot!(F',title="Function and derivative",label="derivative",legend=true)

    plot(F,f,label="function")
    df(x) = -sin(3*x^2)*6*x
    plot!(F',df,label="derivative",legend=true)

    basis = FourierBasis(21,-1,1)⊗FourierBasis(21,-1,1)
    domain = intersect(disk(0.8),disk(0.4))
    f(x,y) = cos(7*x-2*y^2)
    F = Fun(f, basis, domain)

    plot(domain, n=100)

    plot(F,size=(600,400))

    contourf(F,colorbar=true)

    plot(F,plot_ext=true)

    P = contour(F,f)
end

test_plots()

end
