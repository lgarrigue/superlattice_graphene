include("../1_common_graphene_manips/common_graphene.jl")

#################### First step : produce Bloch eigenfunctions of the 3d monolayer, Vks and Vint. The long step is the computation of Vint

function produce_bloch_functions_and_potentials()
    p = get_standard_param_module("dftk_new")
    # Fixed parameters
    p.basis.a = 4.66 # length of the vectors of the fundamental cell, in Bohr
    p.i_state = 4 # u1 will be the i^th eigenmode (i_state=1 is the ground state), u2 the (i+1)^th, u0 the (i-1)^th

    # Changeable monolayers parameters
    # L : periodicity in z (both for mono and bilayer computations)
    # ecut_over_kd2 : ecut/(kd^2)
    # ecut_over_kd2 = 60; p.basis.L = 100
    # ecut_over_kd2 = 41; p.basis.L = 110
    # ecut_over_kd2 = 60; p.basis.L = 110 # <- ideal
    ecut_over_kd2 = 20; p.basis.L = 30
    ecut_over_kd2 = 60; p.basis.L = 110
    p.spec.ecut = ecut_over_kd2*norm_K_cart(p.basis.a)^2
    p.spec.kgrid = [5,5,1] # for computing the Kohn-Sham potential
    p.spec.kgrid = [3,3,1] # for computing the Kohn-Sham potential
    p.spec.tol_scf = 1e-6
    
    px("(ecut/kD^2,L)=(",ecut_over_kd2,",",p.basis.L,") ; kgrid=",p.spec.kgrid)

    # Initialization of some parameters
    init_params(p)

    #################### Computations
    # SCF for the monolayer
    scfres = scf_graphene(p)

    # Computes the Dirac Bloch functions u1, u2 (at the Dirac point)
    get_dirac_eigenmodes(p)
    get_natural_basis_u1_and_u2(p)
    # test_normalizations(p)
    px("Normalization not tested")

    # Computes the Fermi velocity
    # get_fermi_velocity_with_finite_diffs(4,p) # Computing dE/dk with diagonalizations of H(k), should get 0.380 or 0.381
    records_fermi_velocity_and_fixes_gauge_from_p(p)
    # fermi_velocity_from_scalar_products(p) # displays <ϕ_p,(-i∇) ϕ_j> and <ϕ_p,(-i∇) P ϕ_j>, where P is the parity operator
    # double_derivatives_scalar_products(p) # <∂i ϕk,∂j ϕp> (i,j ∈ {1,2}) and <∂z ϕk,∂z ϕp>

    # Symmetry tests
    test_rot_sym(p)
    test_mirror_sym(p)

    # verifications_operator(p)
    # Exports v, u1, u2
    exports_v_u1_u2(p)


    # Plots v, u0, u1, u2
    if false
        px("Makes plots")
        resolution = 30
        n_motifs = 4
        # rapid_plot(u_fc(1,p),p;n_motifs=n_motifs,name="ϕ1",res=resolution,bloch_trsf=false)
        # rapid_plot(-1*u_fc(2,p),p;n_motifs=n_motifs,name="ϕ2",res=resolution,bloch_trsf=false)

        Gu1_fb = G_fb(u_K(1,p),p)
        Tu2_fb = τ(p.u2_fb,shift_K_M(p.basis),p)

        # Diff
        diff_fb = Gu1_fb + Tu2_fb
        diff_dir = ifft(p.basis,p.spec.K_kpt,diff_fb)
        diff_fc = myfft(diff_dir,p.basis.Vol)

        # Gu1
        Gu1_dir = ifft(p.basis,p.spec.K_kpt,Gu1_fb)
        Gu1_fc = myfft(Gu1_dir,p.basis.Vol)

        # Comparision
        px("Test G Φ1 = -Φ2 ",distance(Gu1_fb,-Tu2_fb))
        px("Test G Φ1 = -Φ2 ",norm(diff_fc))
        px("Test G Φ1 = -Φ2 ",norm(diff_dir))

        # pl = Plots.heatmap(intZ(real.(diff_dir),p))
        # display(pl)

        # rapid_plot(diff_fc,p;n_motifs=n_motifs,name="diff",res=resolution,bloch_trsf=false)
        # rapid_plot(Gu1_fc,p;n_motifs=n_motifs,name="GΦ1",res=resolution,bloch_trsf=false)
        # rapid_plot(p.v_monolayer_fc,p;n_motifs=n_motifs,name="v",res=resolution,bloch_trsf=false)
        # rapid_plot(p.non_local_φ1_fc,p;n_motifs=n_motifs,name="non_local_φ1",res=resolution,bloch_trsf=true)
    end
    p
end

p = produce_bloch_functions_and_potentials()
nothing
