using PGFPlotsX, FrameFun, BasisFunctions, DomainSets, BasisFunctions

f = exp
F1 = Fun(f, ChebyshevT(10))
F2 = expansion(F1)
F3 = Fun(f, extensionframe(ChebyshevT(10),0.0..0.5))
F4 = Fun(f, extensionframe(Fourier(10),0.0..0.5))

for F in (F1,F2)
    Plot(F)

    Plot(F;n=300)

    Axis(F,f)
    Axis(F,F)
    Axis(f,F)
end

Plot(F3;plot_extension=true)

Axis(Plot(F4;plot_complex=true)...)

nothing
