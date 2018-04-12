using BasisFunctions
using StaticArrays
using FrameFun
using Domains



function approximation_step_operators_2d(domain, N1, N2, degr1, degr2)
    B1 = BSplineTranslatesBasis(N1,degr1);
    g1 = BasisFunctions.oversampled_grid(B1,2)
    B2 = BSplineTranslatesBasis(N2,degr2);
    g2 = BasisFunctions.oversampled_grid(B2,2)
    B = B1⊗B2;
    g = g1×g2;

    E1 = CirculantOperator(evaluation_matrix(B1[1],g1)[:])*IndexExtensionOperator(span(B1),gridspace(g1),1:2:length(g1))
    E2 = CirculantOperator(evaluation_matrix(B2[1],g2)[:])*IndexExtensionOperator(span(B2),gridspace(g2),1:2:length(g2))

    G1 = CirculantOperator(E1'E1*[1,zeros(length(g1)-1)...]);
    G2 = CirculantOperator(E2'E2*[1,zeros(length(g2)-1)...]);
    G = G1⊗G2
    E = E1⊗E2

    omega_grid = FrameFun.subgrid(g, domain);
    S = evaluation_operator(span(B), omega_grid).operators[end]

    delta_omega_grid = boundary_grid(g, domain);
    SO2B = restriction_operator(omega_grid, delta_omega_grid);

    boundary_support = boundary_support_grid(B, delta_omega_grid, omega_grid);
    SO2S = restriction_operator(omega_grid, boundary_support);

    BEO = boundary_extension_operator(delta_omega_grid, B);

    DG = inv(G);
    Eo = S*E
    DEo = Eo*DG

    A = Eo;
    Z = DEo;

    I = IdentityOperator(dest(A));
    P = I-A*Z';
    omega_grid, A, Z, P, SO2S, BEO, SO2B
end

function approximation_step_operators_3d(domain, N1, N2, N3, degr1, degr2, degr3)
    B1 = BSplineTranslatesBasis(N1,degr1);
    g1 = BasisFunctions.oversampled_grid(B1,2)
    B2 = BSplineTranslatesBasis(N2,degr2);
    g2 = BasisFunctions.oversampled_grid(B2,2)
    B3 = BSplineTranslatesBasis(N3,degr3);
    g3 = BasisFunctions.oversampled_grid(B3,2)
    B = B1⊗B2⊗B3 ;
    g = g1×g2×g3;

    E1 = CirculantOperator(evaluation_matrix(B1[1],g1)[:])*IndexExtensionOperator(span(B1),gridspace(g1),1:2:length(g1))
    E2 = CirculantOperator(evaluation_matrix(B2[1],g2)[:])*IndexExtensionOperator(span(B2),gridspace(g2),1:2:length(g2))
    E3 = CirculantOperator(evaluation_matrix(B3[1],g3)[:])*IndexExtensionOperator(span(B3),gridspace(g3),1:2:length(g3))

    G1 = CirculantOperator(E1'E1*[1,zeros(length(g1)-1)...]);
    G2 = CirculantOperator(E2'E2*[1,zeros(length(g2)-1)...]);
    G3 = CirculantOperator(E3'E3*[1,zeros(length(g3)-1)...]);
    G = G1⊗G2⊗G3
    E = E1⊗E2⊗E3

    omega_grid = FrameFun.subgrid(g, domain);
    S = evaluation_operator(span(B), omega_grid).operators[end]

    delta_omega_grid = boundary_grid(g, domain);
    SO2B = restriction_operator(omega_grid, delta_omega_grid);

    boundary_support = boundary_support_grid(B, delta_omega_grid, omega_grid);
    SO2S = restriction_operator(omega_grid, boundary_support);

    BEO = boundary_extension_operator(delta_omega_grid, B);

    DG = inv(G);
    Eo = S*E
    DEo = Eo*DG

    A = Eo;
    Z = DEo;

    I = IdentityOperator(dest(A));
    P = I-A*Z';
    omega_grid, A, Z, P, SO2S, BEO, SO2B
end


function experiment_2d(f, domain, degr1, degr2)
    Ns = [2^k for k in 10:10]
    for i in 2:4
        e = Symbol(string("e",i))
        @eval $e = zeros(length($Ns))
        for j in 0:3
            t = Symbol(string("t",i,j))
            @eval $t = zeros(length($Ns))
        end
    end

    for (i,N) in enumerate(Ns)
        omega_grid, A, Z, P, SO2S, BEO = approximation_step_operators_2d(domain, N, N, degr1, degr2)

        b = sample(omega_grid, f);
        # r[i] = sum(svd(matrix(P*A))[2] .> epsilon)

        r = 5*N
        # x1, t1[i], a1 = @timed system_solve(A, b, epsilon)
        # e1[i] = norm(A*x1-b)
        x2, t20[i], a2 = @timed projection_solve(A, b, epsilon, r, P, Z)
        e2[i] = norm(A*x2-b)
        x3, t30[i], a3 = @timed selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
        e3[i] = norm(A*x3-b)
        x4, t40[i], a4 = @timed sparse_selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
        e4[i] = norm(A*x4-b)

        for m in 3:4
            x = Symbol(string("x",m))
            t1 = Symbol(string("t",m,1))
            t2 = Symbol(string("t",m,2))
            t3 = Symbol(string("t",m,3))
            method = select_method(m)
            @eval begin
                xi1, $t1[$i], = @timed $method(A, b, epsilon, r, P, BEO, SO2S)
                xi2, $t2[$i], = @timed projection_solve2(A, b, Z, xi1)
                println(size(xi1))
                $x,  $t3[$i], = @timed projection_solve3(xi1, xi2)
            end
        end
    end
    tag = string(maximum(Ns))
    for i in 2:4
        e = Symbol(string("e",i))
        @eval writedlm(string("experiment_data/2e",$i,$tag), $e)
        for j in 0:3
            t = Symbol(string("t",i,j))
            @eval writedlm(string("experiment_data/2t",$i,$j,$tag), $t)
        end
    end
    writedlm("experiment_data/2N"tag, Ns)
end

### Start of experiment
function experiment_3d(f, domain, degr1, degr2, degr3)
    Ns = [10,20,30,40,50,]
    for i in [3,4]
        e = Symbol(string("e",i))
        @eval $e = zeros(length($Ns))
        for j in 0:3
            t = Symbol(string("t",i,j))
            @eval $t = zeros(length($Ns))
        end
    end
    tag = string(maximum(Ns))
    for (i,N) in enumerate(Ns)
        println("N ",N)
        omega_grid, A, Z, P, SO2S, BEO = approximation_step_operators_3d(domain, N, N, N, degr1, degr2, degr3)

        b = sample(omega_grid, f);
        x4, t40[i], a4 = @timed sparse_selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
        e4[i] = norm(A*x4-b)
        println(norm(A*x4-b))
        for ii in 3:4
            e = Symbol(string("e",ii))
            @eval writedlm(string("experiment_data/3e",$ii,$tag), $e)
            for j in 0:3
                t = Symbol(string("t",ii,j))
                @eval writedlm(string("experiment_data/3t",$ii,$j,$tag), $t)
            end
        end
        writedlm("experiment_data/3N"tag, Ns)


        x3, t30[i], a3 = @timed selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
        e3[i] = norm(A*x3-b)
        println(norm(A*x3-b))
        for ii in 3:4
            e = Symbol(string("e",ii))
            @eval writedlm(string("experiment_data/3e",$ii,$tag), $e)
            for j in 0:3
                t = Symbol(string("t",ii,j))
                @eval writedlm(string("experiment_data/3t",$ii,$j,$tag), $t)
            end
        end
        writedlm("experiment_data/3N"tag, Ns)

        # xi1, t41[i], = @timed sparse_selection_projection_solve1(A, b, epsilon, r, P, BEO, SO2S)
        # xi2, t42[i], = @timed projection_solve2(A, b, Z, xi1)
        # x,   t43[i], = @timed projection_solve3(xi1, xi2)


    end

    # for i in 2:4
    #     e = Symbol(string("e",i))
    #     @eval writedlm(string("experiment_data/3e",$i,$tag), $e)
    #     for j in 0:3
    #         t = Symbol(string("t",i,j))
    #         @eval writedlm(string("experiment_data/3t",$i,$j,$tag), $t)
    #     end
    # end
    # writedlm("experiment_data/3N"tag, Ns)
end


function system_solve(A, b, epsilon)
    reshape(LAPACK.gelsy!(matrix(A),b,epsilon)[1],size(src(A))...)
end

function projection_solve1(A, b, epsilon, r, P, E, R)
    W = randn(length(src(A)),r);
    AA = apply_multiple(P*A,W)
    y, rr = LAPACK.gelsy!(AA,P*b,epsilon)
    x1 = W*y;
    x1 = reshape(x1, size(src(A))...)
end

function selection_projection_solve1(A, b, epsilon, r, P, E, R)
    AA = matrix(R*P*A*E)
    y1, rr = LAPACK.gelsy!(AA,R*P*b,epsilon)
    x1 = E*y1
end

function sparse_selection_projection_solve1(A, b, epsilon, r, P, E, R)
    AA = sparse_matrix(P*A*E)
    println("sparsity ", nnz(AA)/prod(size(AA)))
    y1 = AA\(P*b)
    x1 = E*y1
end

function select_method(i)
    i==2 && return projection_solve1
    i==3 && return selection_projection_solve1
    i==4 && return sparse_selection_projection_solve1
end

projection_solve2(A, b, Z, x1) = Z'*(b-A*x1);
projection_solve3(x1, x2) = x1 + x2;

function projection_solve(A, b, epsilon, r, P, Z)
    x1 = projection_solve1(A, b, epsilon, r, P, nothing, nothing)
    x2 = projection_solve2(A, b, Z, x1)
    projection_solve3(x1, x2)
end

function selection_projection_solve(A, b, epsilon, P, Z, E, R)
    x1 = selection_projection_solve1(A, b, epsilon, nothing, P, E, R)
    x2 = projection_solve2(A, b, Z, x1)
    projection_solve3(x1, x2)
end


function sparse_selection_projection_solve(A, b, epsilon, P, Z, E, R)
    x1 = sparse_selection_projection_solve1(A, b, epsilon, nothing, P, E, R)
    x2 = projection_solve2(A, b, Z, x1)
    projection_solve3(x1, x2)
end

f2d = (x,y) -> x*(y-1)^2#*exp(x*y)
center = @SVector [.5,.5]
domain2d = disk(.3,center)
center = @SVector [.5,.5,.5]
domain3d = ball(.3,center)
epsilon=1e-10
N1 = 20; N2 = 20; N3 = 10;
degr1 = 1; degr2 = 2; degr3 = 1;

omega_grid, A, Z, P, SO2S, BEO = approximation_step_operators_2d(domain2d,N1,N2,degr1,degr2)

b = sample(omega_grid, f2d);
r = 5*maximum([N1,N2])

x21, = @timed projection_solve1(A, b, epsilon, r, P, BEO, SO2S)
x31, = @timed selection_projection_solve1(A, b, epsilon, r, P, BEO, SO2S)
x41, = @timed sparse_selection_projection_solve1(A, b, epsilon, r, P, BEO, SO2S)

x22, = @timed projection_solve2(A, b, Z, x21)
x32, = @timed projection_solve2(A, b, Z, x31)
x42, = @timed projection_solve2(A, b, Z, x41)

x2, = @timed projection_solve3(x21, x22)
x3, = @timed projection_solve3(x31, x32)
x4, = @timed projection_solve3(x41, x42)
println(norm(A*x2-b));norm(A*x2-b)
println(norm(A*x3-b));norm(A*x3-b)
println(norm(A*x4-b));norm(A*x4-b)


x1, = @timed system_solve(A, b, epsilon)
x2, = @timed projection_solve(A, b, epsilon, r, P, Z)
x3, = @timed selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
x4, = @timed sparse_selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
println(norm(A*x2-b));norm(A*x2-b)
println(norm(A*x3-b));norm(A*x3-b)
println(norm(A*x4-b));norm(A*x4-b)


f3d = (x,y,z) -> x*(y-1)^2*z#*exp(x*y)
omega_grid, A, Z, P, SO2S, BEO = approximation_step_operators_3d(domain3d,N1,N2,N3,degr1,degr2,degr3)


b = sample(omega_grid, f3d);
r = 3*maximum([N1,N2])^2

x21, = @timed projection_solve1(A, b, epsilon, r, P, BEO, SO2S)
x31, = @timed selection_projection_solve1(A, b, epsilon, r, P, BEO, SO2S)
x41, = @timed sparse_selection_projection_solve1(A, b, epsilon, r, P, BEO, SO2S)

x22, = @timed projection_solve2(A, b, Z, x21)
x32, = @timed projection_solve2(A, b, Z, x31)
x42, = @timed projection_solve2(A, b, Z, x41)

x2, = @timed projection_solve3(x21, x22)
x3, = @timed projection_solve3(x31, x32)
x4, = @timed projection_solve3(x41, x42)
println(norm(A*x2-b));norm(A*x2-b)
println(norm(A*x3-b));norm(A*x3-b)
println(norm(A*x4-b));norm(A*x4-b)


x1, = @timed system_solve(A, b, epsilon)
x2, = @timed projection_solve(A, b, epsilon, r, P, Z)
x3, = @timed selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
x4, = @timed sparse_selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
println(norm(A*x2-b));norm(A*x2-b)
println(norm(A*x3-b));norm(A*x3-b)
println(norm(A*x4-b));norm(A*x4-b)

# experiment_2d(f2d, domain2d, degr1, degr2)
# experiment_3d(f3d, domain3d, degr1, degr2, degr3)


omega_grid, A, Z, P, SO2S, BEO = approximation_step_operators_2d(domain2d,N1,N2,degr1,degr2)

b = sample(omega_grid, f2d);


M = matrix(A-A*Z'*A)
sum(M.!=0)
sum(M.!=0)/prod(size(M))


x = M\(P*b)
norm(M*x-P*b)

N = matrix((A-A*Z'*A)*BEO)
sum(N.!=0)
x = (N\(P*b))
norm((A-A*Z'*A)*BEO*x-P*b)

s = svdvals(matrix(A))
scatter(s,yscale=:log10)
L = matrix(A)

[ sum(L[:,k]) for k in 1:size(L,2)] .>1e-10
sum([ sum(L[:,k]) for k in 1:size(L,2)] .>1e-10)
BEO.subindices


K = matrix(A-Z*A'*A)

norm(K-M)


G = matrix(Z'A)
GG = matrix(A'Z)

norm(s-ss)


I = eye(400)
s = svdvals(G)
ss = svdvals(GG)
sss = svdvals(I-G)
ssss = svdvals(L*(I-G))
sssss = svdvals(L)

plot(s,yscale=:log10)
plot!(ss)
plot!(sss)
plot!(ssss)
plot!(sssss)


tol = 1e-10
u,S,v = svd(L)
amask = S .>= .5;sum(amask)
bmask = tol .< S .<.5;sum(bmask)
cmask = (S .< tol); sum(cmask)
Aa = u[:,amask]*Diagonal(S[amask])*v[:,amask]'
Ab = u[:,bmask]*Diagonal(S[bmask])*v[:,bmask]'
Ac = u[:,cmask]*Diagonal(S[cmask])*v[:,cmask]'

W = Aa
L1 = L-W
plot(svdvals(W),yscale=:log10)
plot!(svdvals(L1),yscale=:log10)

K = matrix(Z)
L2 = K'-pinv(W)

plot(svdvals(pinv(W)),yscale=:log10)
plot!(svdvals(L2),yscale=:log10)


plot(svdvals(matrix(A)))
plot!(svdvals(matrix(Z)))
plot!(svdvals(matrix(DG)))


norm(La+Lb+Lc-L)
matrix(Z)-La
svdvals(matrix(Z)-La)
plot(svdvals(matrix(Z)-La),yscale=:log10)
plot!(svdvals(matrix(Z)),yscale=:log10)
plot!(svdvals(matrix(A)),yscale=:log10)
plot(svdvals(Aa)+1e-34,yscale=:log10)
sum(Aa)
u,S,v = svd(G)
tol = 1e-10
mask1 = (S.>1+tol).| (1-tol.> S .> tol)
mask2 =  S .< tol
mask3 = 1+tol .> S .> 1-tol
spy(log10.(1e-14+abs.(v[:,mask1])),layout=(2,2))
spy!(log10.(1e-14+abs.(v)[:,mask2]),subplot=2)
spy!(log10.(1e-14+abs.(v)[:,mask3]),subplot=3)
spy!(log10.(1e-14+abs.(v)[:,:]),subplot=4)
spy(log10.(1e-14+abs.(L)))
spy(log10.(1e-14+abs.(v)))
L
plotlyjs()
