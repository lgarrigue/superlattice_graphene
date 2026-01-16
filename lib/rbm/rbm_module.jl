include("../spectral_stuff.jl")
include("../gram_schmidt.jl")
using Random # for randn

mutable struct EigenRBM
    V # list of vectors being in the variational space. It's an orthonormal basis
    ℓ
    ν
    dim # dimension of the vectors forming the basis
    generalized # true iff we use generalized eigenmode problems
    P
    Pperp
    function EigenRBM(ℓ,ν;gen=false) # gen = true does not work
        R = new()
        R.V = [] # number of vectors in the basis
        R.ℓ, R.ν = ℓ, ν
        R.dim = -1
        R.generalized = gen
        R.P = "" # orthogonal projector
        R.Pperp = "" # 1-R.P
        R
    end
end

function add_vector(v,R;gram_each=false)
    if !R.generalized
        if length(R.V)==0
            R.dim = length(v)
        end
        @assert length(v) == R.dim
        u_perp = gram_schmidt_one_vector(v,R.V)
        u = gram_each ? u_perp : v
        push!(R.V,u)
    else
        push!(R.V,v)
    end
end

# from an EigenRBM with basis (u1,...uℓ) spanning PH to a EigenRBM with basis of the orthogonal space P⟂H. dim_H is the dimension of the full Hilbert space
function orthogonal_RBM(R) 
    dim_H = R.dim # dimension of the total Hilbert space
    Rperp = EigenRBM(-1,0)
    Pperp = get_Pperp(R)
    dim_Hperp = dim_H - nV(R)
    # px("Produces orthogonal RBM dimH=",dim_H," dimH⟂=",dim_Hperp)
    for i=1:dim_Hperp
        ei = zeros(Float64,dim_H)
        ei[i] = 1
        add_vector(Pperp*ei,Rperp)
    end
    # px("size Rperp ",nV(Rperp))
    Rperp
end

# do a Gram-Schmidt the list of functions l to help orthogonalize

function do_Gram_Schmidt_R(R) 
    R.V = gram_schmidt_family(R.V;force=true)
end

nV(R) = length(R.V) # number of vectors in V

function add_vectors(V,R) # adds all the vectors of V
    for j=1:length(V)
        add_vector(V[j],R)
    end
end

# H is a Hamiltonian action or a matrix, gives the approximated eigenfunctions, in the basis R.V
function rbm_approx_eigenmodes(H,R; full=false,tol=1e-5,stops=true) 
    H_act = H
    if H isa AbstractMatrix
        H_act = X -> H*X
    end
    N = nV(R)
    Hrbm = [R.V[i]⋅(H_act(R.V[j])) for i=1:N, j=1:N]
    B = [R.V[i]⋅R.V[j] for i=1:N, j=1:N]
    # Hrbm = zeros(ComplexF64,(N,N))
    # for i=1:N, j=1:N
        # Hrbm[i,j] = R.V[i]⋅(H*R.V[j])
    # end
    # Hrbm = Hermitian(Hrbm)
    # px("Complex part of Hrbm ",norm(imag.(Hrbm)))
    E, C = R.generalized ? gen_eigenmodes(Hrbm,B; tol=tol) : eigenmodes(Hrbm; tol=tol,stops=stops)
    Ψ = []
    for j=1:N
        φ = sum(C[j][k]*R.V[k] for k=1:N)
        if !R.generalized
            @assert abs(norm(φ)-1)<1e-6
        end
        push!(Ψ,φ)
    end
    if full
        R = B^(-1/2)
        return E, Ψ, R*Hrbm*R
    end
    E, Ψ, C
end

function add_random_vectors(n,R) # adds n random vectors
    for k=1:n
        v = copy(randn(R.dim))
        add_vector(v,R)
    end
end

function init_R_level_μ(R,μ,S)
    for k=0:R.ℓ
        ϕk = X("ϕ",(k,μ),S)
        # if α==3
            # v = copy(randn(R.dim))
            # add_vector(v,R)
        # else
        # if α in [4,5,6] || k==0
        add_vector(ϕk,R)
        # else
            # for j=1:10
                # v = copy(randn(length(ϕk)))
                # add_vector(v,R)
            # end
        # end
    end
end

function init_R_from_first_derivatives(R,S)
    for α=0:R.ν
        init_R_level_μ(R,α,S)
    end
end

function test_R(S,R)
    test_orthonormal_basis_R(R)
    @assert S.p["dim"] > nV(R)+5
    # px("R has ",nV(R)," vectors")
end

### Test

# Tests that R.V is an orthonormal basis

function print_norms_basis(R)
    N = nV(R)
    for i=1:N
        x = norm(R.V[i])
        print("Norm i=",i," is ",x," ")
    end
end

function test_orthonormal_basis_R(R;tol = 1e-8)
    if !R.generalized
        test_orthonormal_basis(R.V;tol = tol)
    end
end

# produces the projector over Span(R.V)
create_projector(R) = [R.V[j][i] for i=1:length(R.V[1]), j=1:length(R.V)]
create_full_projector(R) = sum(R.V[k]*R.V[k]' for k=1:length(R.V))

function get_Pperp(R)
    generate_full_projector(R)
    R.Pperp
end

function get_P(R)
    generate_full_projector(R)
    R.P
end

function generate_full_projector(R)
    if R.P=="" || R.Pperp==""
        P = create_full_projector(R)
        R.P = P
        R.Pperp = I-P
    end
end


function pseudo_inverse_RBM(H,R;E="natural",full=false) # (P(cE - PHP)P)^{-1)
    Es, ϕs = rbm_approx_eigenmodes(H,R)
    pseudo_inverse_from_modes(Es,ϕs;z=E,full=full)
end
