# paramcontext.jl


###############################################################################################
# The Param type groups a number of (immutable) parameters of a Fourier extension
###############################################################################################

immutable FParam{N,T}
	n 		::	NTuple{N,Int}	# Degrees of freedom in each dimension
	m 		::	NTuple{N,Int}	# Number of sampling points in each dimension
	l 		::	NTuple{N,Int}	# Number of sampling points in each dimension on extended domain
	t 		::	NTuple{N,T}		# Extension parameter in each dimension
	Ntot	::	Int				# Total number of degrees of freedom
 	Mtot	::	Int				# Total number of sampling points
 	Ltot	::	Int				# Total number of points on extended domain
end

copy(fp::FParam) = FParam(fp.n, fp.m, fp.l, fp.t, fp.Ntot, fp.Mtot, fp.Ltot)

# Compute everything based on n, m and t. Note that t may be adjusted in order to have integer L.
function FParam{N,T}(n::NTuple{N,Int}, m::NTuple{N,Int}, t::NTuple{N,T})
	l = tuple(Int[floor(t[i]*(m[i]-1)) for i=1:N]...)
	t = tuple(T[l[i]/(m[i]-1) for i=1:N]...)
	FParam(n, m, l, t, prod(n), prod(m), prod(l))
end

# Special case for one dimension: user can discard tuples
FParam{T}(n::Int, m::Int, t::T) = FParam((n,), (m,), (t,))

size(fp::FParam) = fp.n


###############################################################################################
# The FContext type groups properties of a Fourier extension on a certain domain
###############################################################################################

type FContext{N,T,D}
	d			:: D 					# D is the Domain type
	dbox		:: FBox{N,T}			# the domain bounding box
	ebox		:: FBox{N,T}			# the extended domain bounding box
	scratch_N	:: Array{Complex{T},N}	# scratch space of length N
	scratch_M	:: Array{Complex{T},N}	# scratch space of length M
	scratch_L	:: Array{Complex{T},N}	# scratch space of length L
	fftplanN							# out-of-place FFT plan of size N
	ifftplanN							# out-of-place IFFT plan of size N
	fftplanL							# in-place FFT plan of size L
	ifftplanL							# in-place IFFT plan of size L
	gridN		:: Grid{N,T}			# Grid of size N on extended domain
	gridM		:: Grid{N,T}			# Grid of size M on domain
	gridL		:: Grid{N,T}			# Grid of size L on extended domain
	maskN		:: Array{Bool,N}		# Boolean mask of the extended domain of length N
	maskM		:: Array{Bool,N}		# Boolean mask of the domain
	maskL		:: Array{Bool,N}		# Boolean mask of the extended domain of length L
	Q 			:: Int					# Total number of points in the masked domain
end


function FContext{N,T,D}(d::D, fp::FParam{N,T})
	dbox = d.box
	ebox = extend(d.box, fp.t)
	fftplanN = plan_fft(zeros(Complex{T},fp.n...),1:N,FFTW.ESTIMATE|FFTW.MEASURE|FFTW.PATIENT)
	ifftplanN = plan_ifft(zeros(Complex{T},fp.n...),1:N,FFTW.ESTIMATE|FFTW.MEASURE|FFTW.PATIENT)
	fftplanL = plan_fft!(zeros(Complex{T},fp.l...),1:N,FFTW.ESTIMATE|FFTW.MEASURE|FFTW.PATIENT)
	ifftplanL = plan_ifft!(zeros(Complex{T},fp.l...),1:N,FFTW.ESTIMATE|FFTW.MEASURE|FFTW.PATIENT)
	scratch_N = zeros(Complex{T}, fp.n)
	scratch_M = zeros(Complex{T}, fp.m)
	scratch_L = zeros(Complex{T}, fp.l)
	gridN = periodicgrid(ebox, fp.n)
	gridM = Grid(dbox, fp.m)
	gridL = periodicgrid(ebox, fp.l)
	maskN = in(gridN, d)
	maskM = in(gridM, d)
	maskL = in(gridL, d)
	Q = sum(maskM)
	FContext(d, dbox, ebox, scratch_N, scratch_M, scratch_L, fftplanN, ifftplanN, fftplanL, ifftplanL, gridN, gridM, gridL, maskN, maskM, maskL, Q)
end

copy(fc::FContext) = FContext(fc.d, fc.dbox, fc.ebox, copy(fc.scratch_N), copy(fc.scratch_M), copy(fc.scratch_L), fc.fftplanN, fc.ifftplanN, fc.fftplanL, fc.ifftplanL, fc.gridN, fc.gridM, fc.gridL, copy(fc.maskN), copy(fc.maskM), copy(fc.maskL), fc.Q )

# Replace the domain of a context.
function similar(fc::FContext, d::AbstractDomain)
	fc = FContext(d, fc.dbox, fc.ebox, copy(fc.scratch_N), copy(fc.scratch_M), copy(fc.scratch_L), fc.fftplanN, fc.ifftplanN, fc.fftplanL, fc.ifftplanL, fc.gridN, fc.gridM, fc.gridL, copy(fc.maskN), copy(fc.maskM), copy(fc.maskL), fc.Q )
	if fc.d == d
		# nothing left to be done
	else
		recomputemasks!(fc)
	end
	fc
end

# Generate suitable FE parameters, based on a domain, degrees of freedom in each direction, extension parameters in each direction and a sampling factor that determines the amount of oversampling
function fe_paramcontext{N,T}(d::AbstractDomain{N,T}, n::NTuple{N}, t, oversampling)
	m = map(x -> int(round(oversampling*x)), n)
	g = Grid(d.box, m)
	Z = in(Grid(d.box, m), d)
	Q = sum(Z)
	Ntot = prod(n)
	ratio = Q / Ntot
	if ratio < oversampling^N-1
  		# do a heuristic guess of the value of M required to have enough points
    	Mratio = Q / prod(m)
    	m = map(x->int(round(oversampling*(Ntot/Mratio)^(1/N))), n)
	end
	fp = FParam(n, m, t)
	fc = FContext(d, fp)
	fp,fc
end

function fe_paramcontext{N,T}(domain::AbstractDomain{N,T}, dof::Int, extension::T, oversampling::T)
	n = tuple(Int[dof for i=1:N]...)
	t = tuple([extension for i=1:N]...)
	fp,fc = fe_paramcontext(domain, n, t, oversampling)
	fp,fc
end




