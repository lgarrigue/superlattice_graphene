include("../1_common_graphene_manips/common_graphene.jl")
include("../../lib/band_diagram.jl")

function produce_band_diagram(fk=1)
    N = 13
    v_intensity = 2
    p = generate_some_2d_param(N,v_intensity;first_k=fk)
    p.type = "mine_2d"

    # Plots band diagram
    bd = BandDiagram()
    K1 = [-2,1]/3; K2 = [-1,2]/3
    Γ = [0.0,0.0]
    Γ2 = 2*K1-K2
    M = Γ2/2
    Klist = [Γ,K1,M]
    Knames = ["Γ","K_1","M"]
    init_klist_knames(Klist,Knames,bd)
    bd.resolution_bands = 55
    bd.energy_center = 2
    bd.energy_amplitude = 2
    bd.plot_all_bands = false
    bd.n_if_action = p.basis.Nfull
    
    bd.k2lin = kred -> X("kcarts_sq",kred,p.basis.stock)
    bd.n_bands = 50
    # bd.n_bands = 100
    add_lattice_vectors(p.basis.a1, p.basis.a2, bd)
    Hop(kred) = action_for_k(kred,p.v_dir,p.basis;coef_kin=1/fk^2)
    σ = spectrum_on_a_path(Hop,bd;print_progress=true,method="k_to_action")
    name = string("test_2d_N_",p.basis.N,"_fk_",fk)
    plot_band_diagram([σ],"Effpot",bd;post_name=name,numbering_bands=true,extension="pdf")
end

function launch_lobpcg()
    N = 50
    V_intensity = 1
    p = generate_some_2d_param(N,V_intensity)
    diag_monolayer_at_k_mine([1,1],p.v_dir,4,p.basis;tol=1e-4)
end


function launch_lobpcg_from_b()
    println("Start lobpcg test from basis")
    b = Basis(2)
    b.N = 50
    V_intensity = 1
    init_cell_vectors(b) 
    init_cell_infinitesimals(b) 
    v_dir = generate_some_sym_potential_dir(V_intensity,b)
    diag_monolayer_at_k_mine(0*[1,1], v_dir, 4, b;tol=1e-4)
    println("End lobpcg test from basis")
end

# launch_lobpcg_from_b()
# launch_lobpcg()
for fk=1:4
    produce_band_diagram(fk)
end
