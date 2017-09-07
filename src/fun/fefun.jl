function FeFun(f::Function; d=1, T = Float64, Gamma=cartesianproduct(interval(-2.,2.),d), Omega=cartesianproduct(interval(-1,1),d), options...)
    set = tensorproduct(map(x->FourierBasis(0,leftendpoint(x), rightendpoint(x), T), elements(Gamma)))


    fun_optimal_N(f, set, Omega; options...)
end
