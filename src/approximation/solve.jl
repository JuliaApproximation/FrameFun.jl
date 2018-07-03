
# Function with equal functionality, but allocating memory
function az_solve(b, A::DictionaryOperator, Zt::DictionaryOperator, util_operators::DictionaryOperator...;
        cutoff = default_cutoff(A), trunc = truncatedsvd_solve, verbose=false, options...)

    P = plunge_operator(A, Zt)
    x2 = trunc(P*b, P*A, util_operators...; cutoff=cutoff, verbose=verbose, options...)
    x1 = Zt*(b-A*x2)
    x1 + x2
end

# Function with equal functionality, but allocating memory
azs_solve(b, A::DictionaryOperator, Zt::DictionaryOperator, RD::DictionaryOperator, SB::DictionaryOperator;
        trunc = restriction_solve, use_plunge=false, options...) =
    az_solve(b, A, Zt, RD, SB; trunc=trunc, use_plunge=use_plunge, options...)

function az_solve(platform::BasisFunctions.Platform, i, f::Function; R=0, options...)
    a = A(platform, i)
    zt = Zt(platform, i)
    s = sampler(platform, i)
    (R == 0) && (R=estimate_plunge_rank(a))
    az_solve(s*f, a, zt; R=R, options...)
end

function azs_solve(platform::BasisFunctions.Platform, i, f::Function; options...)
    a = A(platform, i)
    zt = Zt(platform, i)
    s = sampler(platform, i)
    rd,sb = spline_util_restriction_operators(platform, i)
    azs_solve(s*f, a, zt, rd', sb; options...)
end
