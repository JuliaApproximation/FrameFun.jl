function FeFun(f::Function; T = Float64, approach=fun_optimal_N, Gamma=interval(-T(2),T(2)), Omega=interval(T(-1),T(1)), options...)
    set = FourierBasis(0, leftendpoint(Gamma), rightendpoint(Gamma), T)
    fun_optimal_N(f, set, Omega; options...)
end

function FeFunNd(f::Function, d::Int; T = Float64, Gamma=cartesianproduct(interval(T(-2),T(2)),d), Omega=cartesianproduct(interval(-T(1),T(1)),d), options...)
    if d==1
        warn("use FEFun indstead")
        return
    end
    set = tensorproduct(map(x->FourierBasis(0,leftendpoint(x), rightendpoint(x), T), elements(Gamma)))
    fun_optimal_N(f, set, Omega; max_logn_coefs =min(12, 8^d),options...)
end
