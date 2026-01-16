include("../lib/band_diagram.jl")
include("../1_graphene_dftk/import_export_funs.jl")
include("../1_graphene_dftk/graphene_dftk.jl")

function graphene()
    tol = 1e-5

    N = 27; Nz = 600; gauge_param = 1
    N = 27; Nz = 160; gauge_param = 1
    N = 20; Nz = 120; gauge_param = 1
    # N = 32; Nz = 864; gauge_param = 1

    # Set parameters and makes imports of eigenfunctions
    p = import_u1_u2_V(N,Nz,gauge_param)

    n_bands = 10

    function k2H(k) # too long to call eigen
        ham = get_ham(k,p)
        H_K_block = ham.blocks[1]
        print("forming matrix")
        @time mat = Hermitian(Array(H_K_block)) # 1 second
        print("full diag")
        @time eigen(mat) # 0.25 second
        mat
    end

    function k2E(k)
        ham = get_ham(k,p)
        println("dftk diag")
        @time data = DFTK.diagonalize_all_kblocks(DFTK.lobpcg_hyper, ham, n_bands + 10; n_conv_check=n_bands, tol=tol) # 0.05 second
        sol = DFTK.select_eigenpairs_all_kblocks(data, 1:n_bands+10)
        Es = sol.λ[1]
        Es
    end

    plots_bands(k2E,p)
    nothing
end

function plots_bands(k2E,p)
    # Plots band diagram
    bd = BandDiagram()
    K1 = [-2,1]/3; K2 = [-1,2]/3
    Γ = [0.0,0.0]; Γ2 = 2*K1-K2; M = Γ2/2
    bd.Klist = [Γ,K1,M]
    bd.Knames = ["Γ","K1","M"]
    bd.resolution_bands = 20
    bd.energy_center = 0
    bd.energy_amplitude = 1
    bd.plot_all_bands = false
    add_lattice_vectors(p.a1,p.a2,bd)

    # bd.k2lin = kred -> get_ham(kred;klin=true)
    bd.n_bands = 10
    # σ = spectrum_on_a_path(k2H,bd;print_progress=true,method="k_to_matrix")
    σ = spectrum_on_a_path(k2E,bd;print_progress=true,method="k_to_E")
    plot_band_diagram([σ],"Effpot",bd;post_name="graphene")
end

graphene()
