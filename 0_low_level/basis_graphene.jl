# using Plots#, LaTeXStrings
include("../lib/solvers/lobpcg.jl")
include("../lib/stocks/stock_module.jl")
using DelimitedFiles, FFTW, CairoMakie
px = println

# The functions are represented in a list. Apply reshape to have them in a matrix or cube or other form


# To store the G_plus_K vectors squares
function create_stock(b)
    # 1) Define the different kinds of recursive objects we will need
    # kinds = ["kcarts_vec","kcarts_sq"]

    # 2) Define default parameters values
    p = Dict()
    p["basis"] = b

    # 3) Define how to compute the quantities
    function compute_kcarts_vec(kred,St::Stock)
        bb = p["basis"]
        t = []
        for ik=1:bb.Nfull
            krd = ik2kred_dim(ik,bb)
            l = length(krd)
            last = length(kred)==2 ? 0 : kred[3]
            krr = [kred[1],kred[2],last][1:l]
            kr = krd + krr
            kc = kred2kcart_dim(kr,bb)
            # println(kc)
            push!(t,kc)
        end
        t
    end

    compute_kcarts_sq(kred,St::Stock) = [norm(k)^2 for k in X("kcarts_vec", kred, p["basis"].stock)]

    # 4) Gather computation functions, same order than kinds list
    kcf = [("kcarts_vec",compute_kcarts_vec), ("kcarts_sq",compute_kcarts_sq)]
    
    # 5) Create stock
    St = Stock(kcf,p)
    St
end

mutable struct Basis
    dim
    dims

    # Discretization numbers
    N # N in each direction
    Nz #
    N2d; N3d
    Nfull

    # Direct lattice
    a
    L # physical length of periodicity for the computations for a monolayer, if 3d
    a1; a2
    lattice_2d
    lattice_3d

    # Direct space
    x_axis_cart
    x_axis_red

    # Fourier lattice
    a1_star; a2_star; a3_star
    kred_x 
    kred_z
    kred_grid_2d
    kred_grid_3d
    kred_grid_dim
    kred2ik_dict_2d
    kcart_grid_2d
    kcart_grid_3d
    kz_cart
    kz_cart_sq
    stock

    # Infinitesimal quantities
    dx; dS; dz; dv

    # Volumes, surfaces, distances
    Vol
    sqrtVol
    sqi # 1/sqrt(cell_area)
    cell_area

    # Dirac point
    K_red_2d
    M_red_2d
    K_red_3d
    K_coords_cart

    # Other
    graphene_lattice_orientation

    function Basis(dim)
        b = new()
        b.dim = dim
        b.a = -1
        b.L = -1
        b.Nz = 1
        b.stock = create_stock(b)
        b
    end
end

function graphene_basis_2d(N;a=3)
    b = Basis(2)
    b.N = N
    init_cell_vectors(b,a=a)
    init_cell_infinitesimals(b)
    b
end

# Initialization of graphene cell structures
# when b.dim and b.a and b.L are init
function init_cell_vectors(b;a=4.66) 
    if b.a == -1 # default values
        # px("Generates default parameters")
        b.a = a # a = 4.66 for true graphene
        b.L = 100
    end
    b.cell_area = sqrt(3)*0.5*b.a^2
    b.sqi = 1/sqrt(b.cell_area)
    b.Vol = b.cell_area*b.L
    a1_unit = [1/2;-sqrt(3)/2]
    a2_unit = [1/2; sqrt(3)/2]
    b.a1, b.a2 = b.a.*(a1_unit, a2_unit)

    # Builds dual cell vectors
    b.a1_star, b.a2_star = a_star_from_a(b.a1, b.a2)

    # Builds lattice
    b.lattice_2d = [b.a1 b.a2]

    b.graphene_lattice_orientation = distance(b.a1,rotM(2π/3)*b.a2)<1e-5 || distance(b.a1,rotM(-2π/3)*b.a2)<1e-5 # true if the angle between a1 and a2 is 2π/3, false if it's π/3
    # px("Lattice orientation ",b.graphene_lattice_orientation)

    # Dirac point K
    b.K_red_2d = -[1/3,1/3]

    K1 = [-2,1]/3; K2 = [-1,2]/3
    Γ = [0.0,0.0]
    Γ2 = 2*K1-K2
    M = Γ2/2

    b.M_red_2d = M
    b.K_red_3d = vcat(b.K_red_2d,[0])
    b.K_coords_cart = kred2kcart_2d_to_3d(b.K_red_2d,b)

    b.a3_star = 2π/b.L
    if b.dim == 3
        b.lattice_3d = [b.lattice_2d    zeros(2,1);
                        zeros(1,2)             b.L]
    end
end

norm_K_cart(a) = 4*π/(3*a) # norm of Dirac's momentum

####################################### Coordinates changes

function a_star_from_a(a1,a2) # vector in cartesian direct space to cartesian Fourier space
    lat = 2π*inv([a1 a2]')
    a1_star = lat[1:2,1]
    a2_star = lat[1:2,2]
    a1_star, a2_star
end

function float2int(x;warning=true,name="") # float to int
    y = floor(Int,x+0.5)
    if abs(x-y)>1e-5 && warning
        px("NOT AN INT in ",name," ",x," ",y)
    end
    y
end

# from α(x) (microcell) builds α(x/n) (macrocell), all in Fourier coefs
# seems to be ok at least when the only non vanishing mode is k=0
function scale_2d(ψ,n,b_microcell,b_macrocell)
    bm, bM = b_microcell, b_macrocell
    β = new_fun_2d(bM) # macro
    kmax_M = maximum(abs.(bM.kred_x))
    kmax_m = maximum(abs.(bm.kred_x))
    for ik_m=1:bm.Nfull
        k = ik2kred_2d(ik_m,bm)
        ik_M = kred2ik_2d(n*k,bM)
        if ik_M == -1
            px("MACRO BASIS IS SMALLER THAN MICRO scaled, problem in scale_2d")
        end
        if ik_M!=-1 && norm(n*k) <= kmax_M/2
            β[ik_M] = ψ[ik_m]
        end
        # β[ik_M] = ψ[ik_m]
    end
    β
end

# V from size M to size N
# V is M×M, frequency-sorted, and we want to get the N×N subset of values of V, where N<M. Cutoff of the Fourier functions (C_n)_{-M ≤ n ≤ M} to (C_n)_{-N ≤ n ≤ N} where N ≤ M (exact bounds would need to be adjusted)
# if N>M, adds zeros
function change_N_matrix_old(V,N) 
    M = size(V,1)
    cf = typeof(V[1,1])
    A = zeros(cf,N,N)
    @assert length(size(V))==2
    kN = myfloat2int.(fftfreq(N)*N;name="kN")
    kM = myfloat2int.(fftfreq(M)*M;name="kM")
    for i=1:N, j=1:N
        kiN, kjN = kN[i], kN[j]
        if (kiN in kM) && (kjN in kM)
            IM = k_inv_1d(kiN,M)
            JM = k_inv_1d(kjN,M)
            A[i,j] = V[IM,JM]
        end
    end
    A
end

change_N_matrix(V,basis_micro,basis_macro) = scale_2d(V,1,basis_micro,basis_macro)

####################################### Coordinates changes

cart2red_mat(M,b) = (b.lattice_2d')*M*inv(b.lattice_2d') # matrix of cartesian Fourier coords to reduced Fourier coords

function cart2red_four(a,c,b) # G = ma^* to m, vector in cartesian Fourier space to reduced Fourier space
    S = cart2red_mat([a c],b)
    (S[1:2,1],S[1:2,2])
end

# (function k_inv_1d) from [k,l] which are momentums in reduced coordinates, ∈ ℤ^2, gives the label [ki,kl] so that C^D_{ki,li} (stored in numerics) = C_{k,l} (true coefficient of computations)
k_inv_1d(k,N) = Int(mod(k,N))+1 # from k in reduced Fourier coordinate to ki such that f^D[ki] = f_k, where f_k = (1/sqrt |Ω|) int e^{-ikx} f(x) dx and f^D is the array storing the coefficients f_k, k = fftfreq[ki] so k_inv_1d inverts fftfreq
k_inv(k,l,b) = (k_inv_1d(k,b.N),k_inv_1d(l,b.N))
kz_inv(k,b) = k_inv_1d(k,b.Nz)
k_inv_v(K,b) = (k_inv_1d(K[1],b.N),k_inv_1d(K[2],b.N))

####################################### Handles graphene cell structure

# Direct map
ik2kred_dim(ik,b) = b.kred_grid_dim[ik]
ik2kred_2d(ik,b) = b.kred_grid_2d[ik]
ik2kred_3d(ik,b) = b.kred_grid_3d[ik]

# Inverse k maps
function init_kred2ik_dict_2d(b)
    b.kred2ik_dict_2d = Dict()
    for ik=1:b.N^2
        kred = b.kred_grid_2d[ik]
        b.kred2ik_dict_2d[kred] = ik
    end
    # px(b.kred2ik_dict_2d)
    # test_kred2ik_2d(b)
end

function kred2ik_2d(kred,b) 
    if kred in b.kred_grid_2d
        return b.kred2ik_dict_2d[kred]
    end
    -1
end

function test_kred2ik_2d(b) 
    for ik=1:b.N^2
        kred = ik2kred_2d(ik,b)
        ikp = kred2ik_2d(kred,b) 
        @assert ik==ikp
    end
end

# Cartesian
ik2kcart_2d(ik,b) = b.kcart_grid_2d[ik]
ik2kcart_3d(ik,b) = b.kcart_grid_3d[ik]

function kred2kcart_2d(kred,b;shift=[0,0])
    k = kred + shift
    k[1]*b.a1_star + k[2]*b.a2_star
end

function kred2kcart_3d(kred,b;shift=[0,0,0])
    k = kred + shift
    k[1]*vcat(b.a1_star,[0]) + k[2]*vcat(b.a2_star,[0]) + k[3]*b.a3_star*[0,0,1]
end

function kred2kcart_dim(kred,b)
    if b.dim==2
        return kred2kcart_2d(kred,b)
    end
    kred2kcart_3d(kred,b)
end

kred2kcart_2d_to_3d(kred_2d,b;shift=[0,0]) = vcat( kred2kcart_2d(kred_2d,b,shift=shift) , [0] )
# kred2kcart_3d_to_3d(kred_3d,b) = kred2kcart_2d(kred_2d,b,shift=shift) , [0] )

# Runs after init_cell_vectors, and one knows the discretization resolution. Initializes other quantities
function init_cell_infinitesimals(b::Basis) # needs a, N, Nz ; only micro quantities !
    # k stuff
    b.kred_x = float2int.(fftfreq(b.N)*b.N;name="kred_x")
    b.kred_z = float2int.(fftfreq(b.Nz)*b.Nz;name="kred_z")
    @assert b.N == length(b.kred_x)
    b.kred_grid_2d = [[b.kred_x[ikx],b.kred_x[iky]] for ikx=1:b.N, iky=1:b.N]
    b.kred_grid_3d = [[b.kred_x[ikx],b.kred_x[iky],b.kred_z[ikz]] for ikx=1:b.N, iky=1:b.N, ikz=1:b.Nz]
    b.kred_grid_dim = b.dim==2 ? b.kred_grid_2d : b.kred_grid_3d
    init_kred2ik_dict_2d(b)

    # Cartesian
    b.kcart_grid_2d = [kred2kcart_2d(b.kred_grid_2d[ikx,iky],b) for ikx=1:b.N, iky=1:b.N]
    b.kcart_grid_3d = [kred2kcart_3d(b.kred_grid_3d[ikx,iky,ikz],b) for ikx=1:b.N, iky=1:b.N, ikz=1:b.Nz]
    b.kz_cart = (2π/b.L) .* b.kred_z
    b.kz_cart_sq = b.kz_cart.^2

    # Other
    b.dims = b.dim==2 ? (b.N,b.N) : (b.N,b.N,b.Nz)
    b.N2d = b.N^2
    b.N3d = b.N2d*b.Nz
    b.Nfull = b.dim == 2 ? b.N2d : b.N3d
    b.dS = b.cell_area/b.N^2 # surface element

    b.dx = b.a/b.N
    b.x_axis_cart = (0:b.N-1)*b.dx
    b.dz = b.L/b.Nz
    b.dv = b.Vol/(b.N^2*b.Nz)
end

kcart_sq(b) = [norm(b.kcart_grid_3d[ikx,iky,ikz])^2 for ikx=1:b.N, iky=1:b.N, ikz=1:b.Nz]

####################################### Fourier transforms

# If a_i = f(x_i) are the true values of function f in direct space, then myfft(f) gives the true Fourier coefficients
# where (𝔽f)_m := 1/sqrt(|Ω|) ∫_Ω f(x) e^{-ima^*x} dx are those coefficients, actually myfft(a)[m] = (𝔽f)_{m-1}
myfft(f,Vol) = FFTW.fft(f)*sqrt(Vol)/length(f)
myifft(f,Vol) = FFTW.ifft(f)*length(f)/sqrt(Vol)

####################################### Scalar products, integrations

function scaprod(ϕ,ψ,b::Basis,four=true) # scalar products
    d = length(size(ϕ))
    @assert d==length(size(ψ))
    if four
        return ϕ⋅ψ
    else
        dVol = d==1 ? b.dx : d==2 ? dS : b.dv
        return dVol*ϕ⋅ψ
    end
end
norm2(ϕ,b,four=true) = real(scaprod(ϕ,ϕ,b,four))
norms(ϕ,b,four=true) = sqrt(norm2(ϕ,b,four))

function Kinetic(u_four,b::Basis) # kinetic energy of u
    (∂1u,∂2u,∂3u) = ∇(u_four,b)
    norm2_3d(∂1u,b)+norm2_3d(∂2u,b)+norm2_3d(∂3u,b)
end

# u in linear, v in Fourier cube, to v*u in linear
function mul(v,u,b;first_is_direct_space=false)
    u_fc = myreshape(u,b)
    if first_is_direct_space # v is in Fourier cube form
        return mylinearize(myfft(v.*myifft(u_fc,b.Vol),b.Vol))
    end
    mylinearize(myfft(myifft(v,b.Vol).*myifft(u_fc,b.Vol),b.Vol))
end

function shift_K_R(b)
    [float2int.((I-matrix_rot_red(2π/3,b))*b.K_red_2d;name="shiftKR");0]
end

function shift_K_M(b)
    M_four_2d = float2int.(cart2red_mat([1 0;0 -1],b);name="M") # mirror symmetry, M_four_2d = [0 1;1 0]
    [float2int.((I-M_four_2d)*b.K_red_2d;name="shiftKM");0]
end

# ############# Universal
# (∇+ik)f where f is in any linear representation (Fourier cube lin or Fourier ball lin)
function ∇_k_low_level(f_lin,G_plus_k_cart,dim) 
    n = length(f_lin)
    Tuple([im*f_lin[i]*G_plus_k_cart[i][a] for i=1:n] for a=1:dim)
end
# (∇+ik)^2 f where f is in any linear representation (Fourier cube lin or Fourier ball lin)
Δ_k(f_lin,G_plus_k_cart_sq) = -f_lin .* G_plus_k_cart_sq

# ############# Mines
∇_k_mine(f_lin,b;kred=[0,0,0]) = ∇_k_low_level(f_lin, X("kcarts_vec",kred,b.stock), b.dim)
Δ_k_mine(f_lin,b;kred=[0,0,0]) = Δ_k(f_lin, X("kcarts_sq",kred,b.stock))

# Laplacian
# -Δf = ∑_{m,m_z} ((ma^*)^2 + (2π/L m_z)^2) hat(f)_m
# the following gives Δ, not -Δ
myreshape(u,b) = reshape(u,b.dims)
mylinearize(u) = vcat(u...)

function Δ_old(f_fc,b;kred=[0,0])
    g = similar(f_fc)
    kx, ky = kred[1], kred[2]
    for m=1:b.N, n=1:b.N, q=1:b.Nz
        kr = b.kred_grid_2d[m,n]
        m0, n0 = kr[1], kr[2]
        k2 = norm((m0+kx)*b.a1_star + (n0+ky)*b.a2_star)^2 + b.kz_cart_sq[q]
        c = f_fc[m,n,q]
        g[m,n,q] = k2*c
    end
    -g
end

# v is the potential, in Fourier cube
apply_schro(u_lin,v_dir_fc,b;kred=[0,0],coef_kin=1) = -0.5*coef_kin*Δ_k_mine(u_lin,b;kred=kred) + mul(v_dir_fc,u_lin,b;first_is_direct_space=true)

kin_energy(u_lin,b;kred=[0,0],coef_kin=1) = -0.5*coef_kin*real(u_lin⋅Δ_k_mine(u_lin,b;kred=kred))

function action_for_k(kred,v_direct_fc,b;coef_kin=1)
    u_lin -> apply_schro(u_lin,v_direct_fc,b;kred=kred,coef_kin=coef_kin)
end
# Gives the proportions of each Fourier mode
function proportions_fourier_modes(ψ,dim,b;err=1e-4,only_kz=false,name="")
    props = Dict()
    n = dim==2 ? b.N2d : b.N3d
    f = dim==2 ? ik2kcart_2d : ik2kcart_3d
    # f = dim==2 ? ik2kred_2d : ik2kred_3d
    for ik=1:n
        k = f(ik,b)
        x = only_kz ? k[3] : k
        nk = round(norm(x),digits=3)
        v = abs(ψ[ik])
        if nk in keys(props)
            props[nk] += v
        else
            props[nk] = v
        end
    end
    norms, prop = [], []
    for k in keys(props)
        push!(norms,k)
        push!(prop,props[k])
    end
    p = sortperm(norms)
    norms = norms[p]
    propr = round.(prop[p],digits=3)
    px("Proportions of Fourier modes : ")
    for j=1:length(norms)
        if sum(prop[j:end])<err
            break
        end
        print("|k|=",norms[j]," prop=",propr[j]," ")
    end
    fig = Figure()
    ax = Axis(fig[1, 1])
    CairoMakie.scatter!(ax,norms, propr, markersize = 3)
    save("props_norms_"*name*".pdf", fig)
    px("\n")
    props
end

function plot_first_modes(ψ,b;cut=3)
    for ik=1:b.Nfull
        kred = ik2kred_dim(ik,b)
        v = ψ[ik]
        # if abs(v-ψ[2])<0.000001
        if norm(kred)<cut
            x = round(abs(v),digits=2)
            print("k=",kred,"v=",x," ")
        end
    end
end

function cut(u,n,b)
    φ = copy(u)
    for ik=1:b.N3d
        kc = ik2kcart_3d(ik,b)
        if norm(kc) > n
            φ[ik] = 0
        end
    end
    φ
end
