# function FeFun(f::Function; T = Float64, approach=fun_optimal_N, Gamma=Interval(-T(2),T(2)), Omega=Interval(T(-1),T(1)), options...)
#     set = Fourier(0, leftendpoint(Gamma), rightendpoint(Gamma), T)
#     fun_optimal_N(f, set, Omega; options...)
# end

function FeFun(f::Function, d::Int=1; T = Float64, Omega=FeFun_Omega_default(T,d), Gamma=2*boundingbox(Omega), options...)
    if d==1
        set = Fourier(0, infimum(Gamma), supremum(Gamma))
    else
        set = tensorproduct(map((x,y)->Fourier(0,x, y),infimum(Gamma),supremum(Gamma)))
    end
    fun_optimal_N(f, set, Omega; max_logn_coefs =min(12, 8^d),options...)
end

FeFun_Omega_default(::Type{T}, d::Int) where {T} = cartesianproduct(Interval(-T(1),T(1)),d)
