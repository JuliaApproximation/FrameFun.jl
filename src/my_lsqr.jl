import IterativeSolvers

realtype{T}(::Type{Complex{T}}) = T
realtype{T}(::Type{T}) = T

import IterativeSolvers.Adivtype
import IterativeSolvers.zerox
import IterativeSolvers.MatrixCFcn
import IterativeSolvers.ConvergenceHistory

import Base: push!

Arealdivtype(A,b) = realtype(typeof(one(eltype(b))/one(eltype(A))))

# A complexified version of the lsqr routine in IterativeSolvers
# which is in turn adapted from the Matlab implementation at
#    http://www.stanford.edu/group/SOL/software/lsqr.html
function my_lsqr!(x, ch::ConvergenceHistory, A, b; damp=0, atol=sqrt(eps(Arealdivtype(A,b))), btol=sqrt(eps(Arealdivtype(A,b))), conlim=one(Arealdivtype(A,b))/sqrt(eps(Arealdivtype(A,b))), maxiter::Int=max(size(A,1), size(A,2)))
    # Sanity-checking
    m = size(A,1)
    n = size(A,2)
    length(x) == n || error("x should be of length ", n)
    length(b) == m || error("b should be of length ", m)
    for i = 1:n
        isfinite(x[i]) || error("Initial guess for x must be finite")
    end

    # Initialize
    empty!(ch)
    ch.threshold = (atol, btol, conlim)
    T = Adivtype(A, b)
    Tr = Arealdivtype(A, b)
    itn = istop = 0
    ctol = conlim > 0 ? convert(Tr,1/conlim) : zero(Tr)
    Anorm = Acond = ddnorm = res2 = xnorm = xxnorm = z = sn2 = zero(T)
    cs2 = -one(T)
    dampsq = damp*damp
    tmpm = Array(T, m)
    tmpn = Array(T, n)

    # Set up the first vectors u and v for the bidiagonalization.
    # These satisfy  beta*u = b-A*x,  alpha*v = A'u.
    u = b-A*x
    v = zeros(T, n)
    beta = norm(u)
    alpha = zero(T)
    if beta > 0
        scale!(u, one(T)/beta)
        Ac_mul_B!(v,A,u)
        alpha = norm(v)
    end
    if abs(alpha) > zero(Tr)
        scale!(v, one(T)/alpha)
    end
    w = copy(v)
    ch.mvps += 2

    Arnorm = alpha*beta
    if Arnorm == 0
        return
    end

    rhobar = alpha
    phibar = bnorm = rnorm = r1norm = r2norm = beta

    #------------------------------------------------------------------
    #     Main iteration loop.
    #------------------------------------------------------------------
    while itn < maxiter
        itn = itn + 1

        # Perform the next step of the bidiagonalization to obtain the
        # next beta, u, alpha, v.  These satisfy the relations
        #      beta*u  =  A*v  - alpha*u,
        #      alpha*v  =  A'*u - beta*v.
        A_mul_B!(tmpm, A, v)
        for i = 1:m
            u[i] = tmpm[i] - alpha*u[i]
        end
        beta = norm(u)
        if beta > 0
            scale!(u, one(T)/beta)
            Anorm = sqrt(Anorm*Anorm + alpha*alpha + beta*beta + dampsq)
            Ac_mul_B!(tmpn, A, u)
            for i = 1:n
                v[i] = tmpn[i] - beta*v[i]
            end
            alpha  = norm(v)
            if alpha > 0
                for i = 1:n v[i] /= alpha; end
            end
        end
        ch.mvps += 2

        # Use a plane rotation to eliminate the damping parameter.
        # This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
        rhobar1 = sqrt(rhobar*rhobar + dampsq)
        cs1     = rhobar/rhobar1
        sn1     = damp  /rhobar1
        psi     = sn1*phibar
        phibar  = cs1*phibar

        # Use a plane rotation to eliminate the subdiagonal element (beta)
        # of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
        rho     =   sqrt(rhobar1*rhobar1 + beta*beta)
        cs      =   rhobar1/rho
        sn      =   beta   /rho
        theta   =   sn*alpha
        rhobar  = - cs*alpha
        phi     =   cs*phibar
        phibar  =   sn*phibar
        tau     =   sn*phi

        # Update x and w
        t1      =   phi  /rho
        t2      = - theta/rho
        for i = 1:n
            wi = w[i]
            x[i] += t1*wi
            w[i] = v[i] + t2*wi
            wirho = wi/rho
            ddnorm += wirho*wirho
        end

        # Use a plane rotation on the right to eliminate the
        # super-diagonal element (theta) of the upper-bidiagonal matrix.
        # Then use the result to estimate  norm(x).
        delta   =   sn2*rho
        gambar  = - cs2*rho
        rhs     =   phi - delta*z
        zbar    =   rhs/gambar
        xnorm   =   sqrt(xxnorm + zbar^2)
        gamma   =   sqrt(gambar*gambar + theta*theta)
        cs2     =   gambar/gamma
        sn2     =   theta /gamma
        z       =   rhs   /gamma
        xxnorm  =   xxnorm + z^2

        # Test for convergence.
        # First, estimate the condition of the matrix  Abar,
        # and the norms of  rbar  and  Abar'rbar.
        Acond   =   Anorm*sqrt(ddnorm)
        res1    =   phibar^2
        res2    =   res2 + psi^2
        rnorm   =   sqrt(res1 + res2)
        Arnorm  =   alpha*abs(tau)

        # 07 Aug 2002:
        # Distinguish between
        #    r1norm = ||b - Ax|| and
        #    r2norm = rnorm in current code
        #           = sqrt(r1norm^2 + damp^2*||x||^2).
        #    Estimate r1norm from
        #    r1norm = sqrt(r2norm^2 - damp^2*||x||^2).
        # Although there is cancellation, it might be accurate enough.
        r1sq    =   rnorm^2 - dampsq*xxnorm
        r1norm  =   sqrt(abs(r1sq));   if real(r1sq) < 0 r1norm = - r1norm; end
        r2norm  =   rnorm
        push!(ch, r1norm)

        # Now use these norms to estimate certain other quantities,
        # some of which will be small near a solution.
        test1   =   abs(rnorm /bnorm)
        test2   =   abs(Arnorm/(Anorm*rnorm))
        test3   =   abs(one(T)/Acond)
        t1      =   abs(test1/(one(T) + Anorm*xnorm/bnorm))
        rtol    =   abs(btol + atol*Anorm*xnorm/bnorm)

        # The following tests guard against extremely small values of
        # atol, btol  or  ctol.  (The user may have set any or all of
        # the parameters  atol, btol, conlim  to 0.)
        # The effect is equivalent to the normal tests using
        # atol = eps,  btol = eps,  conlim = 1/eps.
        if itn >= maxiter  istop = 7; end
        if 1 + test3  <= 1 istop = 6; end
        if 1 + test2  <= 1 istop = 5; end
        if 1 + t1     <= 1 istop = 4; end

        # Allow for tolerances set by the user
        if  test3 <= ctol  istop = 3; end
        if  test2 <= atol  istop = 2; end
        if  test1 <= rtol  istop = 1; end

        if istop > 0 break end
    end
    x
end

function my_lsqr!(x, A, b; kwargs...)
    T = Arealdivtype(A, b)
    z = zero(T)
    ch = ConvergenceHistory(false, (z,z,z), 0, T[])
    my_lsqr!(x, ch, A, b; kwargs...)
    x, ch
end

my_lsqr(A, b; kwargs...) = my_lsqr!(zerox(A, b), A, b; kwargs...)

