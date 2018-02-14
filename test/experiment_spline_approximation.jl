using BasisFunctions
using StaticArrays
using FrameFun
using Domains



function approximation_step_operators(domain, N1, N2, degr1, degr2)
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
    # SO2B = restriction_operator(omega_grid, delta_omega_grid);

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
    omega_grid, A, Z, P, SO2S, BEO
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
    AA = matrix(R*P*A*E)
    AA[abs.(AA).<1e-14] = 0
    SA = sparse(AA)
    y1 = SA\(R*P*b)
    x1 = E*y1
end

function select_method(i)
    i==2 && return projection_solve1
    i==3 && return selection_projection_solve1
    i==4 && return sparse_selection_projection_solve1
end

projection_solve2(A, b, Z, x1) = Z'*(b-A*x1);
projection_solve3(x1, x2) = x1 + x2;

function projection_solve(A, b, epsilon, R, P, Z)
    W = randn(length(src(A)),R+5);
    AA = apply_multiple(P*A,W)
    y, rr = LAPACK.gelsy!(AA,P*b,epsilon)
    x1 = W*y;
    x1 = reshape(x1, size(src(A))...)
    x2 = Z'*(b-A*x1);
    x1+x2;
end

function selection_projection_solve(A, b, epsilon, P, Z, E, R)
    AA = matrix(R*P*A*E)
    y1, rr = LAPACK.gelsy!(AA,R*P*b,epsilon)
    x1 = E*y1
    x2 = Z'*(b-A*x1);
    x1+x2;
end

function sparse_selection_projection_solve(A, b, epsilon, P, Z, E, R)
    AA = matrix(R*P*A*E)
    AA[abs.(AA).<1e-14] = 0
    SA = sparse(AA)
    y1 = SA\(R*P*b)

    x1 = E*y1
    x2 = Z'*(b-A*x1);
    x1+x2;
end

f = (x,y) -> x*(y-1)^2*exp(x*y)
center = @SVector [.5,.5]
domain  = disk(.3,center)

epsilon=1e-10
N1 = 40; N2 = 40;
degr1 = 1; degr2 = 2;

omega_grid, A, Z, P, SO2S, BEO = approximation_step_operators(domain,N1,N2,degr1,degr2)

b = sample(omega_grid, f);
r = sum(svd(matrix(P*A))[2] .> epsilon)

x21, = @timed projection_solve1(A, b, epsilon, r, P, BEO, SO2S)
x31, = @timed selection_projection_solve1(A, b, epsilon, r, P, BEO, SO2S)
x41, = @timed sparse_selection_projection_solve1(A, b, epsilon, r, P, BEO, SO2S)

x22, = @timed projection_solve2(A, b, Z, x21)
x32, = @timed projection_solve2(A, b, Z, x31)
x42, = @timed projection_solve2(A, b, Z, x41)

x2, = @timed projection_solve3(x21, x22)
norm(A*x2-b)
x3, = @timed projection_solve3(x31, x32)
norm(A*x3-b)
x4, = @timed projection_solve3(x41, x42)
norm(A*x4-b)




x1, = @timed system_solve(A, b, epsilon)
norm(A*x1-b)
x2, = @timed projection_solve(A, b, epsilon, r, P, Z)
norm(A*x2-b)
x3, = @timed selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
norm(A*x3-b)
x4, = @timed sparse_selection_projection_solve(A, b, epsilon, P, Z, BEO, SO2S)
norm(A*x4-b)

### Start of experiment
Ns = [2^k for k in 3:9]
for i in 2:4
    e = Symbol(string("e",i))
    @eval $e = zeros(length(Ns))
    for j in 0:3
        t = Symbol(string("t",i,j))
        @eval $t = zeros(length(Ns))
    end
end

for (i,N) in enumerate(Ns)
    omega_grid, A, Z, P, SO2S, BEO = approximation_step_operators(domain, N, N, degr1, degr2)

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

    for m in 2:4
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

Nsextended = [Ns...,(100:100:2000)...]
# step 2,3 complexity
tstep2 = zeros(length(Nsextended))

tstep3 = zeros(length(Nsextended))
for (i,N) in enumerate(Nsextended)
    omega_grid, A, Z, P, SO2S, BEO = approximation_step_operators(domain, N, N, degr1, degr2)
    b = sample(omega_grid, f);
    x1 = rand(N,N)
    x2 = rand(N,N)
    xx,tstep2[i] = @timed projection_solve2(A, b, Z, x1)
    @assert size(xx) == size(x1)
    _,tstep3[i] = @timed projection_solve3(x1, x2)
end

tag = string(maximum(Ns))
for i in 2:4
    e = Symbol(string("e",i))
    @eval writedlm(string("experiment_data/e",$i,$tag), $e)
    for j in 0:3
        t = Symbol(string("t",i,j))
        @eval writedlm(string("experiment_data/t",$i,$j,$tag), $t)
    end

writedlm("experiment_data/N"tag, Ns)
writedlm("experiment_data/Ne"tag, Nsextended)
writedlm("experiment_data/tstep2",tstep2)
writedlm("experiment_data/tstep3",tstep3)
