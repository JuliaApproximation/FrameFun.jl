#module ExpFun

import Base
using ArrayViews
import IterativeSolvers



# Generate suitable FE parameters, based on a domain, degrees of freedom in each direction, extension parameters in each direction and a sampling factor that determines the amount of oversampling
function create_FE_setting(d::Domain, N1=10, N2=10, T1=2, T2=2, samplingfactor=2)
	M1 = round(samplingfactor*(2*N1+1))
	M2 = round(samplingfactor*(2*N2+1))
	Z = evaluate_on_grid(d, Grid(M1,M2))
	Q = sum(Z)
	Ntot = (2*N1+1)*(2*N2+1)
	ratio = Q / Ntot
	if ratio < samplingfactor^2-1
  		# do a heuristic guess of the value of M required to have enough points
    	Mratio = Q / (M1*M2)
    	M1 = round(samplingfactor*sqrt(Ntot/Mratio))
    	M2 = M1
	end
  
	L1 = floor(T1 * (M1-1))
	L2 = floor(T2 * (M2-1))
  
	fp = FE_Param(M1, M2, L1, L2, N1, N2)
	fc = FE_Context(d, fp)
	fp,fc
end


function FE_matrix(fp::FE_Param, fc::FE_Context)
  N1 = fp.N1; N2 = fp.N2; L1 = fp.L1; L2 = fp.L2; M1 = fp.M1; M2 = fp.M2; N = fp.N;
  
  A = zeros(Complex128, fc.Q, N)
  r = zeros(Complex128, N)
  u = zeros(Complex128, fc.Q)

  r[1] = 1.0
  FE_matrixvectorproduct!(r, u, fp, fc)
  A[:,1] = u
  
  for i=2:N
    r[i-1] = 0.0
    r[i] = 1.0
    FE_matrixvectorproduct!(r, u, fp, fc)
    A[:,i] = u
  end
  return A
end

function reshape_N_to_L!(c, d, N1, N2, L1, L2)

	# First zero out d
	d[1:end] = 0
	d[1:N1+1,1:N2+1] = view(c, 1:N1+1, 1:N2+1);
	d[1:N1+1,end-N2+1:end] = view(c, 1:N1+1, 2*N2+1-N2+1:2*N2+1);
	d[end-N1+1:end,1:N2+1] = view(c, 2*N1+1-N1+1:2*N1+1,1:N2+1);
	d[end-N1+1:end,end-N2+1:end] = view(c, 2*N1+1-N1+1:2*N1+1, 2*N2+1-N2+1:2*N2+1);
end

function select_entries_from_mask!(y,z,mask)
	l = 1
	for j=1:size(mask,2)
		for i=1:size(mask,1)
			if mask[i,j]
				y[l] = z[i,j]
				l+=1
			end
		end
	end
end

# Do the matrix-vector product with an FE matrix and place the output in y
function FE_matrixvectorproduct!(x, y, fp::FE_Param, fc::FE_Context)

	N1 = fp.N1; N2 = fp.N2; L1 = fp.L1; L2 = fp.L2; M1 = fp.M1; M2 = fp.M2; N = fp.N;
  
	c = reshape(x, 2*N1+1, 2*N2+1)
	d = fc.scratch_L
	reshape_N_to_L!(c, d, N1, N2, L1, L2)
	fc.fftplanL(d)
  
	z = fc.scratch_M
#	z[:] = view(d, 1:M1, 1:M2)
	for j=1:M2
		for i=1:M1
			z[i,j] = d[i,j]
		end
	end
#	y[:] = z[fc.maskM]
	select_entries_from_mask!(y,z,fc.maskM)
end


function FE_rhs(f, fp::FE_Param, fc::FE_Context)
  N1 = fp.N1; N2 = fp.N2; L1 = fp.L1; L2 = fp.L2; M1 = fp.M1; M2 = fp.M2; N = fp.N;
  
  (x,y) = create_grid_points(fc.box, fp)

  z = zeros(Complex128, M1, M2)
  for i=1:length(x)
    for j=1:length(y)
      z[i,j] = f(x[i], y[j])
    end
  end
  B = z[fc.maskM]
end

function pseudo_backslash(A,B)
    F = qrfact(A, pivot=true)
    Q = full(F[:Q],thin=true)
    R = F[:R]
    p = F[:p]
    Rd = abs(diag(R))
    threshold = 2e-11
    K = minimum(find(Rd .< threshold))-1
    R1 = R[:,1:K]
    x1 = R1 \ (Q'*B)
    x = zeros(Complex128,size(A,2))
    x[p[1:K]] = x1
    x
end


function truncated_svd_solve(A,B,threshold=2e-11)
	(u,s,v) = svd(A)
	I = minimum(find(s .< threshold))
	s = s[1:I]
	u = u[:,1:I]
	v = v[:,1:I]
	y = diagm(s) \ (u'*B)
	x = v*y
end


function FE_approximate(f, fp, fc)
	A = FE_matrix(fp, fc)
	B = FE_rhs(f, fp, fc)
	q,r = qr(A)
	c = r \ (q'*B)
#	c = truncated_svd_solve(A, B)
#	c = A \ B
	(c,A,B)
end


function FE_approximate_iterative(f, fp, fc)
  A = FE_matrix(fp, fc)
  B = FE_rhs(f, fp, fc)
  (c,ch) = IterativeSolvers.lsqr(A,B)
  (c,A,B,ch)
end

function relocate_point(x, y, b::Box)
        u = (x .- left(b))/(right(b)-left(b))
        v = (y .- bottom(b))/(top(b)-bottom(b))
        u,v
end


function FE_evaluate(x, y, c, fp::FE_Param, fc::FE_Context)
	N1 = fp.N1; N2 = fp.N2; L1 = fp.L1; L2 = fp.L2; M1 = fp.M1; M2 = fp.M2; N = fp.N;
	T1 = fp.T1; T2 = fp.T2;

	x,y = relocate_point(x,y,fc.box)
	d = reshape(c, 2*N1+1, 2*N2+1)
	freq_i = [0:N1, -N1:-1]'
	freq_j = [0:N2, -N2:-1]'

	z = zeros(Complex128, size(x))
	for i=1:length(x)
		F = exp(-2*pi*1im*freq_i'*x[i]/T1) * exp(-2*pi*1im*freq_j*y[i]/T2)
		z[i] = sum(d .* F)
	end
	z
end

FE_evaluate(x::Real, y::Real, c, fp::FE_Param, fc::FE_Context) = (z=FE_evaluate([x], [y], c, fp, fc); z[1])




#end
