include("../1_common_graphene_manips/common_graphene.jl")

include("../../lib/band_diagram.jl")
include("../../lib/pt/asymptotic_behavior.jl")
include("../../lib/rbm/test_rbm_pt.jl")

include("../2_monolayer_derivatives/graphene_derivatives.jl")
include("../3_effective_system/stock_module_coefs.jl")
include("../3_effective_system/eff_op.jl")
include("../3_effective_system/solve_exact_scaled.jl")
include("../3_effective_system/generate_parameters.jl")

#################### Second step : compute the effective potentials 𝕍, 𝕎, 𝔸, etc. Very rapid

function test()
    N = 15
    λ = 0.1
    # Basis
    b = Basis(2)
    b.N = N
    init_cell_vectors(b) 
    init_cell_infinitesimals(b)
    ############## -Δ+v
    # Parameters
    p = Params(2)
    p.basis = b

    # Potential
    v_intensity = 2
    p.v_dir = generate_some_sym_potential_dir(v_intensity, b; fourier=false, first_k=1)

    # Fills eigenmodes at K
    get_dirac_eigenmodes_mine(p)
    px("Energies ",p.Es_K[1:8])

    #Other
    p.i_state = 1
    p.type = "mine_2d"

    # Rotate eigenmodes in the right way
    get_u1_u2_in_good_basis_from_p(p)


    px("\n\n-Δ+v+λV")
    ############## -Δ+v+λV
    # Parameters
    p2 = Params(2)
    p2.basis = b

    # Potential
    V = generate_some_sym_potential_dir(1,b)
    p2.v_dir = p.v_dir + λ*V

    # Fills eigenmodes at K
    get_dirac_eigenmodes_mine(p2)
    px("Energies ",p2.Es_K[1:8])

    #Other
    p2.i_state = 1
    p2.type = "mine_2d"

    # Rotate eigenmodes in the right way
    get_u1_u2_in_good_basis_from_p(p2)
    


    # Exact
    ϕ_ex = p2.us_K[p2.i_state+2]

    # Approx
    ϕ0 = p.us_K[p.i_state+2]
    Vϕ0 = mul_3d(V,ϕ0,b)
    ϕ1 = apply_pseudo_inverse(Vϕ0,p)
    ϕ_approx = ϕ0+λ*ϕ1

    di = dist_quo(ϕ_ex,ϕ_approx)
    px("di ",di)
end

test()
