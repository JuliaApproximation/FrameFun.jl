# fourierextension.jl
# This file groups functions related to the construction and evaluation of Fourier extensions.

function frequencyvector(n) 
	nh = int((n-1)/2)
	[0:nh, -nh:-1]
end


# Embed the matrix c in an otherwise zero matrix d (overwrites d)
@ngenerate N nothing function pad_with_zeros!{N,T}(c::Array{T,N}, d::Array{T,N})
	fill!(d,zero(T))

	@nloops N i c begin
		(@nref N d i) = @nref N c i
	end
end

# Truncate a larger matrix d to a smaller matrix c (overwrites c)
@ngenerate N nothing function truncate!{N,T}(c::Array{T,N}, d::Array{T,N})
	@nloops N i c begin
		(@nref N c i) = @nref N d i
	end
end


# Reshape functions: we want to efficiently copy the data from a vector of length N to a larger vector of length L.
# Hard to do with Cartesian, but it can be done for any dimension recursively:
function reshape_N_to_L!{N}(c, d, n::NTuple{N}, l::NTuple{N})
	nh = map(x->int((x-1)/2), n)
	# First zero out d
	fill!(d, 0)
	reshape_N_to_L_rec!(c, d, (), (), nh, n, l)
	reshape_N_to_L_rec!(c, d, (), (), nh, n, l)
end

function reshape_N_to_L_rec!{N}(c, d, c_ranges, d_ranges, nh::NTuple{N}, n, l)
	reshape_N_to_L_rec!(c, d, tuple(c_ranges...,1:nh[1]+1), tuple(d_ranges...,1:nh[1]+1), nh[2:end], n[2:end], l[2:end])
	reshape_N_to_L_rec!(c, d, tuple(c_ranges...,n[1]-nh[1]+1:n[1]), tuple(d_ranges...,l[1]-nh[1]+1:l[1]), nh[2:end], n[2:end], l[2:end])
end

# The end of the recursion: perform the actual copy
function reshape_N_to_L_rec!(c, d, c_ranges, d_ranges, nh::NTuple{1}, n::NTuple{1}, l::NTuple{1})
	# Currently, the two lines below do some allocation. Using views is not a great improvement.
	# d[d_ranges...,1:nh[1]+1] = c[c_ranges...,1:nh[1]+1]
	# d[d_ranges...,l[1]-nh[1]+1:l[1]] = c[c_ranges...,n[1]-nh[1]+1:n[1]]
	copy_ranges!(c, d, tuple(c_ranges...,1:nh[1]+1), tuple(d_ranges...,1:nh[1]+1))
	copy_ranges!(c, d, tuple(c_ranges...,n[1]-nh[1]+1:n[1]), tuple(d_ranges...,l[1]-nh[1]+1:l[1]))
end


function reshape_L_to_N!{N}(c, d, n::NTuple{N}, l::NTuple{N})
	nh = map(x->int((x-1)/2), n)
	reshape_L_to_N_rec!(c, d, (), (), nh, n, l)
	reshape_L_to_N_rec!(c, d, (), (), nh, n, l)
end

function reshape_L_to_N_rec!{N}(c, d, c_ranges, d_ranges, nh::NTuple{N}, n, l)
	reshape_L_to_N_rec!(c, d, tuple(c_ranges...,1:nh[1]+1), tuple(d_ranges...,1:nh[1]+1), nh[2:end], n[2:end], l[2:end])
	reshape_L_to_N_rec!(c, d, tuple(c_ranges...,n[1]-nh[1]+1:n[1]), tuple(d_ranges...,l[1]-nh[1]+1:l[1]), nh[2:end], n[2:end], l[2:end])
end

# The end of the recursion: perform the actual copy
function reshape_L_to_N_rec!(c, d, c_ranges, d_ranges, nh::NTuple{1}, n::NTuple{1}, l::NTuple{1})
	copy_ranges!(d, c, tuple(d_ranges...,1:nh[1]+1), tuple(c_ranges...,1:nh[1]+1))
	copy_ranges!(d, c, tuple(d_ranges...,l[1]-nh[1]+1:l[1]), tuple(c_ranges...,n[1]-nh[1]+1:n[1]))
end

# Perform the copy without additional allocation
@ngenerate N nothing function copy_ranges!{N}(c, d, c_ranges::NTuple{N}, d_ranges::NTuple{N})
	@nloops N i x->1:length(c_ranges[x]) begin
		(@nref N d x->d_ranges[x][i_x]) = (@nref N c x->c_ranges[x][i_x])
	end
end



# Implements y = z[mask]
@ngenerate N nothing function maskselect!{N}(y, z, mask::Array{Bool,N})
	l = 1
	@nloops N i mask begin
		if (@nref N mask i)
			y[l] = (@nref N z i)
			l += 1
		end
	end
end

# Implements z[mask] = y
@ngenerate N nothing function maskfill!{N}(y, z, mask::Array{Bool,N})
	l = 1
	@nloops N i mask begin
		if (@nref N mask i)
			(@nref N z i) = y[l]
			l += 1
		end
	end
end


# Compute y = Ax. Overwrites y.
# x has dimensions of Fourier coefficients (n...).
# y has dimensions of length of the mask
@ngenerate N nothing function matrixvectorproduct!{N,T}(x, y, fp::FParam{N,T}, fc::FContext{N,T})
	n = fp.n;
	m = fp.m;
	t = fp.t;
	l = fp.l;
  
	d = fc.scratch_L
	c = reshape(x, n...)
	reshape_N_to_L!(c, d, n, l)
	fc.fftplanL(d)

	# Copy below is not really needed when using maskselect, it is when using y[:]=z[mask]   
	# z = fc.scratch_M
	# @nloops N i z begin
	# 	(@nref N z i) = (@nref N d i)
	# end
	# y[:] = d[fc.maskM]
	maskselect!(y,d,fc.maskM)
end

# Compute x = A'y. Overwrites x.
@ngenerate N nothing function matrixvectorproduct_transpose!{N,T}(x, y, fp::FParam{N,T}, fc::FContext{N,T})
	n = fp.n;
	m = fp.m;
	t = fp.t;
	l = fp.l;
  
	d = fc.scratch_L
	fill!(d, 0)
	maskfill!(y, d, fc.maskM)
	fc.ifftplanL(d)

	c = reshape(x, n...)
	reshape_L_to_N!(c, d, n, l)
	x[:] = x*prod(l)
end


function fe_matrix{N,T}(fp::FParam{N,T}, fc::FContext{N,T})
	Ntot = fp.Ntot;
	Mtot = fp.Mtot;
 	A = zeros(Complex{T}, fc.Q, Ntot)
 	r = zeros(Complex{T}, Ntot)
 	u = zeros(Complex{T}, fc.Q)

 	r[1] = one(T)
 	matrixvectorproduct!(r, u, fp, fc)
 	A[:,1] = u

 	for i=2:Ntot
		r[i-1] = zero(T)
 		r[i] = one(T)
 		matrixvectorproduct!(r, u, fp, fc)
		A[:,i] = u
	end
	A
end


@ngenerate N Array{Complex{T},1} function fe_rhs{N,T}(f, fp::FParam{N,T}, fc::FContext{N,T})
	n = fp.n
	m = fp.m
	l = fp.l
	grid = range(fc.gridM);
	z = zeros(Complex{T}, m...)
	@nloops N i x->1:m[x] begin
		(@nref N z i) = f( (@ntuple N x->grid[x][i_x])...)
	end
	z[fc.maskM]
end


function pseudo_backslash(A, B, threshold = 2e-11)
    F = qrfact(A, pivot=true)
    Q = full(F[:Q],thin=true)
    R = F[:R]
    p = F[:p]
    Rd = abs(diag(R))
    I = find(Rd .< threshold)
    if length(I) > 0
    	K = minimum(I)-1
    else
    	K = size(R,2)
    end
    R1 = R[:,1:K]
    x1 = R1 \ (Q'*B)
    x = zeros(eltype(A),size(A,2))
    x[p[1:K]] = x1
    x
end


function fe_approximate_direct(f, fp, fc)
	A = fe_matrix(fp, fc)
	B = fe_rhs(f, fp, fc)
	c = pseudo_backslash(A, B)
	c = reshape(c, fp.n...)
	(c,A,B)
end

function fe_approximate_lsqr{N,T}(f, fp::FParam{N,T}, fc)

	# Define functions appropriate for the IterativeSolvers framework
	my_A_mul_B!(output, x) = (matrixvectorproduct!(x, output, fp, fc); output)

	my_Ac_mul_B!(output, y) = (matrixvectorproduct_transpose!(output, y, fp, fc); output)

	B = fe_rhs(f, fp, fc)
	y,ch = my_lsqr(MatrixCFcn{Complex{T}}(fc.Q, fp.Ntot, my_A_mul_B!, my_Ac_mul_B!), B, maxiter = 100)
	c = reshape(y, fp.n...)
	print("Stopped after ", ch.mvps, " iterations with residual ", abs(ch.residuals[end]), ".")
	c
end



# function truncated_svd_solve(A,B,threshold=2e-11)
# 	(u,s,v) = svd(A)
# 	I = minimum(find(s .< threshold))
# 	s = s[1:I]
# 	u = u[:,1:I]
# 	v = v[:,1:I]
# 	y = diagm(s) \ (u'*B)
# 	x = v*y
# end



# function FE_approximate_iterative(f, fp, fc)
#   A = FE_matrix(fp, fc)
#   B = FE_rhs(f, fp, fc)
#   (c,ch) = IterativeSolvers.lsqr(A,B)
#   (c,A,B,ch)
# end

@ngenerate N Complex{T} function evaluate{N,T}(fp::FParam{N,T}, fc::FContext{N,T}, coef, x::AbstractVector)
	n = fp.n;
	t = fp.t;
	nh = map( x-> int((x-1)/2), n)

	@nexprs N j->freq_j = [0:nh[j], -nh[j]:-1]
	z = 0
	@nloops N i j->1:length(freq_j) begin
		z = z + (@nref N coef i) * prod( @ntuple N j->exp(-2*pi*1im*freq_j[i_j]*x[j]/t[j]))
	end
	z
end


# Evaluate on a grid. The grid has to be on the extended domain!
function evaluate!{N,T}(fp::FParam{N,T}, fc::FContext{N,T}, coef, l::NTuple{N}, result)
	n = fp.n;
	reshape_N_to_L!(coef, result, n, l)
	fft!(result)
end

function evaluate{N,T}(fp::FParam{N,T}, fc::FContext{N,T}, coef, l::NTuple{N})
	z = zeros(Complex{T}, l...)
	evaluate!(fp, fc, coef, l, z)
	z
end


function differentiate{N,T}(c::AbstractArray{T,N}, j)
	n = size(c,j)
	freq = frequencyvector(n)
	d = c .* reshape(freq, ones(j-1)..., n, ones(N-j))
end



