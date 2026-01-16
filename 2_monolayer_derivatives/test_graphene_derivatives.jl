include("../1_graphene_dftk/graphene_dftk.jl")
include("../1_graphene_dftk/import_export_funs.jl")
include("../1_common_graphene_manips/common_graphene.jl")
# include("../lib/rbm/rbm_module.jl")
include("../../lib/rbm/test_rbm_pt.jl")
include("../../lib/pt/asymptotic_behavior.jl")
include("../../lib/linalg/linear_independence.jl")

include("../2_monolayer_derivatives/graphene_derivatives.jl")
include("./manip_tensors.jl")

# tests if vectors are linearly independent
function test_independence()
    type = "dftk_import"
    type = "mine_2d"
    p = get_standard_param_module(type;N=15)
    ℓ = 2
    de = Derivatives(ℓ,p;n_zeroth=3,di_triplet=2)
    l = all_states(de)

    px("Number of states : ",length(l))
    test_linear_independence(l;tol=1e-3)
end

function compute_matrices()
    type = "dftk_import"
    # type = "mine_2d"
    p = get_standard_param_module(type;N=27)
    # p = get_standard_param_module(type;N=15)
    test_rot_sym(p)
    test_M_sym(p)
    ℓ = 2
    de = Derivatives(ℓ,p;n_zeroth=2,di_triplet=2)
    print_first_order(de)
    print_second_order(de)
end

function test_derivative()
    # type = "dftk_import"
    type = "mine_2d"
    dk = [1,-2,0.0]
    k_dep = "k_dependent"
    p = get_standard_param_module(type;N=9)
    test_PC_sym(p)
    test_M_sym(p)
    test_rot_sym(p)
    px("<ϕ1,(-i∇)ϕ3>=", [u_K(1,p)⋅mi∇(u_K(3,p),p)[j] for j=1:2])
    px("<ϕ3,(-i∇)ϕ3>=", [u_K(3,p)⋅mi∇(u_K(3,p),p)[j] for j=1:2])
    px("<ϕ1,(-i∇)ϕ3> in eV =", [hartree_to_ev*u_K(1,p)⋅mi∇(u_K(3,p),p)[j] for j=1:2])
    

    ℓ = 2
    de = Derivatives(ℓ,p; n_zeroth=2, di_triplet=2)
    plot_Es_around(p)

    # Builds reduced basis stuff
    # is it good that eigenvectors are not orthonormalized ?
    R = EigenRBM(ℓ,0) # builds RBM basis
    # if nV(R) ≤ 17

    states = []
    if k_dep == "k_independent"
        states = all_states(de)
    else
        states = all_directional_states(dk[1:2],de)
    end
    add_vectors(states,R)
    # end
    # do_Gram_Schmidt_R(R;n=20)
    do_Gram_Schmidt_R(R)
    print_norms_basis(R)
    test_orthonormal_basis_R(R;tol = 1e-5)
    px("Number of vectors in basis ",nV(R))

    # Builds λs
    λs = get_λs(-4,0;res=6) # prepares some λ
    # λs = get_λs_loc(-0.4,-0.5,res=5) # prepares some λ
    eigenvecs = []
    eigenvals = []
    # Computes solutions
    for λ in λs
        k = p.basis.K_red_3d + λ*dk
        # exact solution
        action = ""; Es = ""; us = ""
        if type == "dftk_import"
            (Es, us, _, _, mat_at_K, _) = diag_monolayer_at_k_dftk(k,p)
            action = mat_at_K
        else
            n_eigen = 3
            Es, us = diag_monolayer_at_k_mine(k,p.v_dir,n_eigen,p.basis)
            action = action_for_k(k,p.v_dir,p.basis)
        end
        E_rbm, ϕ_rbm, C_rbm = rbm_approx_eigenmodes(action,R) # effective
        E1_exact, E2_exact, u1_exact, u2_exact = extract_eigenmodes(Es,us,p) # exact

        # rbm solution
        u1_rbm, u2_rbm = ϕ_rbm[1], ϕ_rbm[2]
        E1_rbm, E2_rbm = E_rbm[1], E_rbm[2]
        C1, C2 = C_rbm[1], C_rbm[2]
        props = print_props(C1)
        px("dk=",dk)
        if de.ℓ ≥ 1
            px("λ*ε*dk = ",λ*norm(dk)," ratio ",λ*norm(dk)/props[2]) # close to 0.387 = vF !
        end
        # comparision
        # di = dist_quo(φ,u_exact)
        P2(u,v) = u*u' + v*v'
        P_exact = P2(u1_exact, u2_exact)
        P_rbm = P2(u1_rbm, u2_rbm)
        di = dist(P_exact, P_rbm)
        # di = dist(E1_rbm,E1_exact)
        push!(eigenvecs,di)
        dif = abs(E1_rbm - E1_exact)
        push!(eigenvals,dif)
        if dif < 1e-30
            px("Precision insuffisante")
        end
        # props = 
    end
    px("With ",length(states)," states")
    print_results(λs,eigenvecs,1, ℓ+1)
    print_results(λs,eigenvals,1, 2*(ℓ+1))
end

# test_independence()
# compute_matrices()
# println("###############")
test_derivative()
