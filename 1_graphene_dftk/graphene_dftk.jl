# include("stock_module_coefs.jl")
using DFTK, LinearAlgebra, JLD, LaTeXStrings
using Plots, CairoMakie#, DataFrames
setup_threading()
px = println

######################### Parameters

# Class containing all the parameters and functions which need to be stored
mutable struct InfoDFTK
    # Dirac point quantities
    K_kpt # Dirac point in DFTK format

    H # operator H-E
    ham
    H_matrix
    H_K_block
    
    ecut; scfres; dftk_basis
    kgrid
    Gvectors; Gvectors_inv; Gplusk_vectors_cart # Gplusk is used to plot the e^{iKx} u(x) for instance
    # recip_lattice; recip_lattice_inv # from k in reduced coords to a k in cartesian ones
    tol_scf
    function InfoDFTK()
        indf = new()
        indf
    end
end

u_dir(u_fb,p) = DFTK.ifft(p.spec.dftk_basis,p.spec.K_kpt,u_fb)

function u_dir_n(n::Int,p) # n is 0 1 or 2
    u_fb = Dict(0 => u_K(0,p), 1 => u_K(1,p), 2 => u_K(2,p))
    u = u_fb[n]
    u_dir(u,p)
end

function u_fb2fc_dftk(u_fb,p)
    ud = u_dir(u_fb,p)
    myfft(ud,p.basis.Vol)
end

u_fcl(u_fb,p) = mylinearize(u_fb2fc_dftk(u_fb,p))

u_fcl_n(n,p) = mylinearize(myfft(u_dir_n(n,p),p.basis.Vol))
pot_fc(p) = myfft(p.v_dir,p.basis.Vol)

function init_params(p)
    init_cell_vectors(p.basis)

    # Builds some Fourier operations
    R_four_2d = matrix_rot_red(-2π/3,p.basis)

    dftk_Vol = DFTK.compute_unit_cell_volume(p.basis.lattice_3d)
    @assert abs(dftk_Vol-p.basis.Vol) < 1e-10
    # p.Vol_recip = (2π)^3/p.Vol volume of reciprocal lattice
    # p.recip_lattice = DFTK.compute_recip_lattice(p.basis.lattice_3d)
    # p.recip_lattice_inv = inv(p.recip_lattice)

    # Dirac point
    # (1-R_{2π/3})K = m a^*, shift of this vector

    p
end

# function HmE(ξ,p)
    # px("attention H() pas précis. Et HmE lent")
    # -0.5*Δ(ξ, p; kred=p.basis.K_red_2d) + mul_3d(pot_fc(p), ξ, p) - E_fermi(p)*ξ
# end

# px("attention H() pas précis. Et HmE lent")

∇_K_dftk(f_fbl,p) = ∇_k_low_level(f_fbl,p.spec.Gplusk_vectors_cart,3)
Δ_K_dftk(f_fbl,p) = Δ_k(f_fbl, [norm(k)^2 for k in p.spec.Gplusk_vectors_cart] )
H_fb_dftk(ξ,p) = p.spec.H_matrix * ξ

######################### Solves SCF or Schrödinger

# Call SCF to get the Kohn-Sham potential of the monolayer
function scf_graphene(p)
    # Loads model
    psp = load_psp("hgh/pbe/c-q4")
    # Choose carbon atoms
    C = ElementPsp(:C, psp=psp)
    # Puts carbon atoms at the right positions for graphene
    minus = p.basis.graphene_lattice_orientation ? -1 : 1
    C1 = [1/3,minus*1/3,0.0]; C2 = -C1
    # px("C1",C1)
    shifts_atoms = [C1,C2]
    # Builds model
    model = model_PBE(p.basis.lattice_3d, [C,C], shifts_atoms; temperature=1e-3, smearing=Smearing.Gaussian())
    # px("DFTK recip lattice ",model.recip_lattice)
    # px("My recip lattice ",p.a1_star," ",p.a2_star," ",2π/p.basis.L)
    # Builds basis
    dftk_basis = DFTK.PlaneWaveBasis(model; Ecut=p.spec.ecut, p.spec.kgrid)
    # px("Number of spin components ",model.n_spin_components)
    # Loads sizes
    (p.basis.N,p.basis.N,p.basis.Nz) = dftk_basis.fft_size
    # Builds axis
    p.basis.N3d = p.basis.N^2*p.basis.Nz
    # p.x_axis_red = ((0:p.basis.N -1))/p.basis.N
    # p.z_axis_red = ((0:p.basis.Nz-1))/p.basis.Nz
    # p.z_axis_cart = p.z_axis_red*p.basis.L
    # Do initializations
    init_cell_infinitesimals(p.basis)
    @assert abs(p.basis.dv - dftk_basis.dvol) < 1e-10
    # Extracts kpoint information
    i_kpt = 1
    kpt = dftk_basis.kpoints[i_kpt]
    # Runs SCF
    px("Runs DFTK's SCF algorithm for the monolayer...")
    p.spec.scfres = self_consistent_field(dftk_basis)
    px("Energies")
    display(p.spec.scfres.energies)
    px("DFTK's SCF algorithm ended")
    # Extracts the Kohn-Sham potential
    p.v_dir = DFTK.total_local_potential(p.spec.scfres.ham)[:,:,:,1]
    # substract_by_far_value(p.v_dir,p) # NE PAS FAIRE CA !!!

    occupation = floor(Int,sum(p.spec.scfres.ρ)*p.basis.dv +0.5)
    px("Number of effective particles in the model: ",occupation)
    # Extracts Kohn-Sham orbitals
    us = []
    # for i=1:occupation
        # ui = p.spec.scfres.ψ[i_kpt][:, i]
        # push!(us,ui)
        # ui_dir = ifft(dftk_basis,kpt,ui)
    # end
    u0_fb = p.spec.scfres.ψ[i_kpt][:, 1]
    px("Size of initial Fourier ball: ",length(u0_fb))
    px("Size of Fourier cube: ",p.basis.N,"×",p.basis.N,"×",p.basis.Nz)
end


size_b(dftk_basis) = length(r_vectors(dftk_basis))

function basis_from_k(k,p)
    kgrid = ExplicitKpoints([k], [1.])
    dftk_basis = PlaneWaveBasis(p.spec.scfres.basis, kgrid)
    # px("new size ",size_b(dftk_basis)," prev ",size_b(p.spec.scfres.basis))
    dftk_basis
end

function ham_to_mat(ham)
    H_k_block = ham.blocks[1]
    H = Hermitian(Array(H_k_block))
    # px("size ham ",size(H))
    H
end

function get_ham(k,p;klin=false)
    kpt = [k[1],k[2],0.0]
    kk = DFTK.ExplicitKpoints([kpt], [1.])
    # kk = DFTK.ExplicitKpoints([k[1],k[2],0])
    nbasis = DFTK.PlaneWaveBasis(p.spec.scfres.basis, kk)
    K_kpt = nbasis.kpoints[1]
    if klin
        gpk = collect(DFTK.Gplusk_vectors_cart(nbasis, K_kpt))
        gpk = norm.(gpk).^2
        return gpk
    end
    print("forming hamiltonian")
    @time ham = DFTK.Hamiltonian(nbasis; p.spec.scfres.ρ) # 0.01 second
    ham
end

# Compute the band structure at momentum k, when the KS potential is already known
function diag_monolayer_at_k_dftk(k,p;n_bands=1000+p.i_state) # k is in reduced coordinates, n_bands is the number of eigenvalues we want
    # DFTK 0.5.10 version
    # kpt = Kpoint(p.spec.scfres.basis, k, 1)
    # ksymops = [[DFTK.one(DFTK.SymOp)] for _ in 1:length(K_dirac_coord)]
    # ksymops = [[DFTK.identity_symop()] for _ in 1:length(K_dirac_coord)]
    # dftk_basis = PlaneWaveBasis(p.spec.scfres.basis, [kpt])
    # dftk_basis = PlaneWaveBasis(p.spec.scfres.basis, [k], [1])
    
    # DFTK new version
    dftk_basis = basis_from_k(k,p)
    ham = Hamiltonian(dftk_basis; p.spec.scfres.ρ)
    Es, us = "", ""
    ####################### EXACT
    exact_solver = false
    mat = ""
    if exact_solver
        px("Exact solver in graphene_dftk.jl")
        H_matrix = ham_to_mat(ham)
        # px("Hamiltonian size ",size(H_matrix,1))
        (Es,Us) = eigen(H_matrix)
        us = [Us[:,i] for i=1:size(H_matrix)[1]]
        mat = H_matrix
    else
        px("Approximate solver in graphene_dftk.jl")
        data = diagonalize_all_kblocks(DFTK.lobpcg_hyper, ham, n_bands + 3;
                                       n_conv_check=n_bands, tol=p.spec.tol_scf)#, show_progress=true)
        px("Approximate solver done")
        if !data.converged
            @warn "Eigensolver not converged" iterations=data.iterations
        end
        # Extracts solutions
        sol = DFTK.select_eigenpairs_all_kblocks(data, 1:n_bands)
        Es = sol.λ[1]
        us = [sol.X[1][:,i] for i=1:size(sol.X[1],2)]
    end
    ####################### LOBPCG
    # px("Energies are ",Es)

    # Fixes the gauges
    kpt = dftk_basis.kpoints[1]
    ref_vec = ones(p.basis.N,p.basis.N,p.basis.Nz)
    ref_gauge = ref_vec/norm(ref_vec) 
    prv = us[1]
    n_eigenmodes = length(Es)
    for i=1:n_eigenmodes
        ui = us[i]
        ui_dir = ifft(dftk_basis,kpt,ui)
        λ = sum(conj(ui_dir).*ref_gauge)
        c = λ/abs(λ)
        ui_dir *= c
        us[i] = fft(dftk_basis,kpt,ui_dir)
    end
    prv1 = ifft(dftk_basis,kpt,prv)
    prv2 = fft(dftk_basis,kpt,prv1)
    @assert distance(prv,prv2) < 1e-13 # tests if fft ∘ ifft does not loose information

    (Es, us, dftk_basis, ham, mat, kpt)
end

function eigenvecs_at_k(k_red,p)
    (Es_K,us_fb,dftk_basis,ham_kpt) = diag_monolayer_at_k_dftk(k_red,p;n_bands=500)
    u1 = us_fb[p.i_state]
    u2 = us_fb[p.i_state+1]
    u1, u2
end

# Obtains u0, u1 and u2 at the Dirac point
function get_dirac_eigenmodes(p;check_K_cart=true)
    # Do the diagonalization
    (Es_K,us_fb,dftk_basis,ham,kpt) = diag_monolayer_at_k_dftk(p.basis.K_red_3d,p)
    # Extracts basis
    p.us_K = us_fb
    p.Es_K = Es_K
    p.spec.dftk_basis = dftk_basis
    p.spec.K_kpt = p.spec.dftk_basis.kpoints[1]
    K_norm_cart = 4π/(3*p.basis.a)
    if check_K_cart
        k_red2cart_3d(k_red,p) = k_red[1]*vcat(p.basis.a1_star,[0]) + k_red[2]*vcat(p.basis.a2_star,[0]) + k_red[3]*p.basis.a3_star*[0,0,1]
        x = abs(K_norm_cart - norm(k_red2cart_3d(p.spec.K_kpt.coordinate,p)))
        @assert x < 1e-10
    end
    # Matrix and operatros
    px("Computes ham")
    p.spec.ham = ham
    px("Computes HKblock")
    p.spec.H_K_block = ham.blocks[1]
    # px("Computes H_matrix")
    # p.spec.H_matrix = Hermitian(Array(p.spec.H_K_block))
    px("Size of new Fourier ball (at K) ",N_basis(p)," which is the one taken")

    # Extracts dual vectors
    p.spec.Gvectors = collect(G_vectors(p.spec.dftk_basis, p.spec.K_kpt))
    p.spec.Gvectors_inv = inverse_dict_from_array(Tuple.(p.spec.Gvectors)) # builds the inverse
    p.spec.Gplusk_vectors_cart = collect(Gplusk_vectors_cart(p.spec.dftk_basis, p.spec.K_kpt))
end

######################### Coordinates changes

# coords in reduced to coords in cartesian
# k_red2cart_3d(k_red,p) = p.recip_lattice*k_red
# coords in cartesian to coords in reduced
# k_cart2red(k_cart,p) = p.recip_lattice_inv*k_cart

######################### Computes the Fermi velocity

function form_∇_term(u,w,j,p) # <u,(-i∂_j +K_j)w>
    GpKj = [p.spec.Gplusk_vectors_cart[iG][j] for iG=1:N_basis(p)]
    sum((conj.(u)) .* (GpKj.*w))
end

function form_∇2_term(u,w,l,j,p) # <(-i∂_ℓ +K_ℓ)u,(-i∂_j +K_j)w>
    GpKl = [p.spec.Gplusk_vectors_cart[iG][l] for iG=1:N_basis(p)]
    GpKj = [p.spec.Gplusk_vectors_cart[iG][j] for iG=1:N_basis(p)]
    sum((conj.(GpKl.*u)) .* (GpKj.*w))
end

# 2×2 matrices <u,(-i∂_j +K_j)u>
form_∇_one_matrix(u1,u2,j,p) = Hermitian([form_∇_term(u1,u1,j,p) form_∇_term(u1,u2,j,p);
                                          form_∇_term(u2,u1,j,p) form_∇_term(u2,u2,j,p)])

function fermi_velocity_from_scalar_products(p) # Computes the Fermi velocity from the formula <ϕ_p,(-i∇)ϕ_j> = vF σ. Needs that u1 and u2 are rotated before calling
    (A1,A2) = (form_∇_one_matrix(u_K(1,p),u_K(2,p),1,p),form_∇_one_matrix(u_K(1,p),u_K(2,p),2,p))
    # display(A1); display(A2)
    px("Fermi velocity from rotated u1 and u2 ",abs(A1[1,2]))
end

function fermi_velocity_from_scalar_products(p) # <ϕ_i,(-i∇) P ϕ_j> scalar products with parity operation
    u1 = u_K(1,p); u2 = u_K(2,p)
    Pu1 = P_fb(u1,p); Pu2 = P_fb(u2,p)
    v12 = form_∇_term(u1,u2,1,p); v12_ = form_∇_term(u1,u2,2,p)
    vP12 = form_∇_term(u1,Pu2,1,p); vP12_ = form_∇_term(u1,Pu2,2,p)
    vP21 = form_∇_term(u2,Pu1,1,p); vP21_ = form_∇_term(u2,Pu1,2,p)
    c = hartree_to_ev
    px("Fermi velocity <ϕ1,(-i∇) ϕ2> =(",c*v12,",",c*v12_,") eV \n <ϕ1,(-i∇) P ϕ2>",c*vP12," ",c*vP12_," eV, \n<ϕ2,(-i∇) P ϕ1> ",c*vP21," ",c*vP21_," eV")
end

function get_fermi_velocity_with_finite_diffs(n_samplings_over_2,p)
    @assert n_samplings_over_2 ≥ 4
    n_samplings = 2*n_samplings_over_2+1 # has to be even
    Dλ = 0.0001
    start_λ = 1-Dλ; end_λ = 1+Dλ; dλ = (end_λ-start_λ)/(n_samplings-1)
    set_coefs = [start_λ + i*dλ for i=0:n_samplings-1]
    set_cart_K = [λ*p.basis.K_coords_cart for λ in set_coefs]
    values_down = []; values_up = []
    for i=1:n_samplings
        k = k_cart2red(set_cart_K[i],p)
        Es = (diag_monolayer_at_k_dftk(k,p))[1]
        push!(values_down,Es[p.i_state])
        push!(values_up,Es[p.i_state+1])
    end
    Ipoint = Int(floor(n_samplings/4)-1)
    dK = norm(set_cart_K[Ipoint+1]-set_cart_K[Ipoint])
    x_axis = norm.(set_cart_K)
    pl = Plots.plot(x_axis,[values_down,values_up])
    savefig(pl,p.path_plots*"graph_fermi_velocity.png")
    vF = (values_down[Ipoint+1]-values_down[Ipoint])/dK
    vF_up = (values_up[Ipoint+1]-values_up[Ipoint])/dK
    px("Fermi velocity: ",vF,". Verification with upper eigenvalue: ",vF_up)
    vF
end

double_der_matrix(u,w,p) = [form_∇2_term(u,w,i,j,p) for i=1:2, j=1:2]

function double_derivatives_scalar_products(p)
    u1 = u_K(1,p); u2 = u_K(2,p)
    M11 = double_der_matrix(u1,u1,p)
    M12 = double_der_matrix(u1,u2,p)
    c = hartree_to_ev
    px("Double derivatives scalar products in eV\n","<∂i ϕ1,∂j ϕ1> :")
    display(c*M11)
    px("<∂i ϕ1,∂j ϕ2> :")
    display(c*M12)
    px("<∂z ϕ1,∂z ϕ1> : ",c*form_∇2_term(u1,u1,3,3,p))
    px("<∂z ϕ1,∂z ϕ2> : ",c*form_∇2_term(u1,u2,3,3,p))
end



function compute_scaprod_mi∇(u,v,p) # <u,(-i∇_x)v>
    # (∂1_v,∂2_v,∂3_v) = ∇_K_dftk(v,p.basis)
    (∂1_v,∂2_v,∂3_v) = ∇_k_mine(v,p.basis)
    c1 = -im*scaprod(u,∂1_v,p.basis,true); c2 = -im*scaprod(u,∂2_v,p.basis,true)
    c1, c2
end

compute_scaprod_Dirac_u1_u2(p) = compute_scaprod_mi∇(u_fcl_n(1,p),u_fcl_n(2,p),p)

function compute_scaprod_Dirac_u1_u1(p)
    c1,c2 = compute_scaprod_mi∇(u_fcl_n(1,p),u_fcl_n(1,p),p)
    A = kred2kcart_2d(p.basis.K_red_2d,p)
    (c1+A[1],c2+A[2])
end




######################### Test symmetries

function test_normalizations(p)
    # Tests normalization
    # px("Normalization of u1: ",norms(u_dir(1,p),p,false))
    # px("Orthonormality |<u1,u2>|= ",abs(scaprod(u_fcl_n(1,p),u_fcl_n(2,p),p)) + abs(scaprod(u_dir(1,p),p.u2_dir,p,false)))
    @assert abs(norms(u_dir(1,p),p,false)-1) + abs(scaprod(u_fcl_n(1,p),u_fcl_n(2,p),p.basis)) + abs(scaprod(u_dir(1,p),p.u2_dir,p.basis,false)) < 1e-5
end


function test_mirror_sym(p)
    px("Tests for mirror symmetry, Gf(x1,x2,z) := f(x1,-x2,z)")
    (Gu0,Tu0) = (M_fb(u_K(0,p),p),τ_fb(u_K(0,p),shift_K_M(p.basis),p))
    (Gu1,Tu1) = (M_fb(u_K(1,p),p),τ_fb(u_K(1,p),shift_K_M(p.basis),p))
    (Gu2,Tu2) = (M_fb(u_K(2,p),p),τ_fb(u_K(2,p),shift_K_M(p.basis),p))

    px("Test Φ1(-z)=-Φ1(z) ",distance(parity_z(u_fcl_n(1,p),p.basis),-u_fcl_n(1,p)))
    px("Test Φ2(-z)=-Φ2(z) ",distance(parity_z(u_fcl_n(2,p),p.basis),-u_fcl_n(2,p)))

    px("Test Φ1(-z)=Φ1(z) ",distance(parity_z(u_fcl_n(1,p),p.basis),u_fcl_n(1,p)))
    px("Test Φ2(-z)=Φ2(z) ",distance(parity_z(u_fcl_n(2,p),p.basis),u_fcl_n(2,p)))

    px("Test G Φ0 =  Φ0 ",distance(Gu0,Tu0))
    px("Test G Φ0 = -Φ0 ",distance(Gu0,-Tu0))

    px("Test G Φ0 = -Φ1 ",distance(Gu0,-Tu1))
    px("Test G Φ0 =  Φ1 ",distance(Gu0,Tu1))

    px("Test G Φ0 = -Φ2 ",distance(Gu0,-Tu2))
    px("Test G Φ0 =  Φ2 ",distance(Gu0,Tu2))

    px("Test G Φ1 = -Φ2 ",distance(Gu1,-Tu2))
    px("Test G Φ1 =  Φ2 ",distance(Gu1,Tu2))

    px("Test G v  =   v ",distance(pot_fc(p),M_four(pot_fc(p),p.basis)))
end
