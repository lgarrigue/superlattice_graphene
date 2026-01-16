include("../../lib/plot/plot_matrix.jl")
COEF_v_tot = 1

function modulo_K(n,b)
    k = b.K_red_2d
    if n % 3 == 0
        k = [0,0]
    elseif n % 3 == 2
        k = -b.K_red_2d
    end
    k
    # mod(n,3)*b.K_red_2d
end


# number of eigenmodes
function N_eigenmodes_exact(invε,p_micro)
    # map = Dict(6=>60, 7=>75, 8=>95, 9=>116, 10=> 140, 11=>165, 12=>196, 13=>225, 14=>270, 15=>290, 16=>335, 17=>380, 18=>420, 19=>465, 20=>515)
    P(x) = 1.222*x^2 + 0.77*x + 18 # interpolation polynomial
    # if invε==7
        # return 80
    # elseif invε==8
        # return 95
    # end
    if invε<=8
        # ep = ifo.no_V ? 1 : 1/invε
        ep = 1/invε
        m = floor(Int,(1/ep)^2*(p_micro.i_state+3ep+1)+2)
        # if ifo.no_V
            # return m + 5
        # end
        return m
    else
        return floor(Int,P(invε))
    end
end

# size of basis
function N_basis_exact(invε,N_micro_basis;precision="normal")
    # if ifo.no_V
        # return N_micro
    # end
    # px("~~~~~~~ ",invε," ",N_micro)
    fact = 1
    if precision == "precise"
        fact = 2
    elseif precision == "no"
        fact = 0.5
        if invε ≤ 5 fact = 1 end # otherwise the basis is too small
    end
    n = floor(Int,N_micro_basis*invε*fact)
    px("n basis exact ",n)
    n
end

function compute_exact_potential(coords,St_sol)
    invε, precision, V_intensity, V_kind = coords
    # Extract quantities
    be = X("basis_exact",(invε,precision),St_sol)
    p = X("p_micro",precision,St_sol)
    N_macro_local = prec2N(precision)
    bm = X("basis_macro",N_macro_local,St_sol)
    V_dir = X("V_dir",(V_intensity,V_kind,N_macro_local),St_sol)
    # Microscopic potential
    v_fc = myfft(p.v_dir, p.basis.Vol)
    v_scaled_fc = scale_2d(v_fc, invε, p.basis, be)
    v_scaled_dir = myifft(v_scaled_fc, be.Vol)
    # Macroscopic potential
    V_fc = myfft(V_dir,bm.Vol)
    V_new_basis = change_N_matrix(V_fc, bm, be)
    V_dir = myifft(V_new_basis, be.Vol)
    # Total potential
    v_tot = v_scaled_dir + (1/invε)*V_dir
    if false # Check that v and V did not change their scales with the manipulations
        Mv = maximum(abs.(v_scaled_dir))
        MV = maximum(abs.(V_dir))
        px("v/V ",Mv/MV)
    end
    # plot_mat(v_scaled_dir;name="pot",dir="plots_wavefuns")
    v_f = myfft(v_scaled_dir, be.Vol)
    # plot_mat(v_f;name="pot_four",dir="plots_wavefuns")
    # print("v plotted ")
    # vflin = mylinearize(v_f)
    # test_nΩ_periodicity(vflin,invε,be)
    v_tot
end

function get_v_exact_dir(oc,St_sol)
    coords = oc.invε, oc.precision, oc.V_intensity, oc.V_kind
    X("v_ex_dir",coords,St_sol)
end

function solve_exact(oc, St_sol; n_eigenmodes_min=5, tol=1e-6, maxiter=300, V_equal_0=false, k_exact=true, verbose=true)
    t0 = time()
    be = get_basis_exact(oc,St_sol)
    p = get_micro_p(oc,St_sol)
    # k_exact : if takes k/ε (exact one) or k/ε mod 3 (approximates). ε=1/18 : 1mn15 en k_exact et 42sc en modulo. Be careful with this when one wants eigenfunctions !
    ε = 1/oc.invε
    Es, ψs = "", ""
    p = get_micro_p(oc,St_sol)
    # if p.type[1:4]=="dftk" 
    # px("CETTE BRANCHE EST OBSOLETE")
    # K3d = vcat(k_red_shifted, [0])
    # (Es, us_exact, basis_kred, _, _, kpt) = diag_monolayer_at_k_dftk(K3d,p)
    # ψ_lin = ψs[p.i_state]
    # ψ_dir = DFTK.ifft(basis_kred,kpt,ψ_lin)
    # ψ_ex = myfft(ψ_dir,p.basis.Vol)
    # else
    # coef_kin = ifo.no_V ? 1 : ε^2 # Kinetic coefficient
    k_red_shifted = ""

    # if ifo.no_V
    if false
        k_red_shifted = p.basis.K_red_2d + ε*oc.kred
    else
        # if k_exact
        k_red_shifted = p.basis.K_red_2d*oc.invε + oc.kred
            # k_red_shifted = mod(oc.invε,3)*p.basis.K_red_2d*oc.invε + oc.kred
        # else
            # k_red_shifted = modulo_K(oc.invε,p.basis) + oc.kred # k
        # end
    end
    n_eigenmodes = max(n_eigenmodes_min, N_eigenmodes_exact(oc.invε,p)) # Number of eigenmodes
    FULL_DIAG = oc.invε ≤ 4 ? true : false # ε=1/18 : 2mn20 en full et 1mn15 en lobpcg
    v_tot = COEF_v_tot*get_v_exact_dir(oc,St_sol)
    Es, ψs = diag_monolayer_at_k_mine(k_red_shifted, v_tot, n_eigenmodes, be ;
                                      tol=tol, maxiter=maxiter, coef_kin=ε^2, full_diag=FULL_DIAG,prints=MAKE_TESTS)

    # KEEP THIS
    # if !k_exact  # !ifo.no_V && 
        # t = p.basis.K_red_2d*oc.invε + oc.kred - k_red_shifted 
        # for i=1:length(ψs)
            # ψs[i] = τ_fc_2d(ψs[i],-t,be)
            # ψs[i] /= norm(ψs[i])
        # end
    # end
    # get_natural_basis_u1_and_u2(p)
    # test_R_inv(myfft(p.v_dir,p.basis.N),p.basis;name="v_dir")
    # test_rot_sym(p)
    # records_fermi_velocity_and_fixes_gauge(p) 
    # end
    Es = real.(Es)
    Es = ε^(-1)*(Es .- E_fermi(p))
    check_at_least_two_positive(Es)
    if verbose
        print("Solved exact in ",round(time()-t0,digits=3)," sec. ")
    end
    Es, ψs
end

function plot_ψs_ex(k0,Es0,ψs,oc,St_sol)
    ε = 1/oc.invε
    be = get_basis_exact(oc,St_sol)
    i = get_i_dmode(oc,St_sol)
    if norm(k0)≤1e-6
        px("Energies ",Es0)
        for j=-10:10
            I = i+j
            if 1 ≤ I ≤ length(ψs)
                ψ = ψs[I]
                # φ = myreshape(ψ, be)
                u = abs.(FFTW.ifft(ψ))
                # px(typeof(ψ)," len ",length(ψ)," ",length(X("kcarts_sq",k0,be.stock)))
                # lφ = mylinearize(ψ)
                # kin = kin_energy(lφ,be;kred=k0,coef_kin=ε^2)
                # px("Kin of ",I,"=",kin)
                name = string(oc.invε,"_",I,"_te_",round(Es0[I],digits=3),"_kin_",round(0,digits=3))
                plot_mat(u;name=name,dir="plots_wavefuns")
                u_four = abs.(ψ)
                plot_mat(u_four;name=name*"_four",dir="plots_wavefuns")
                px("i_dirac is ",i)
                # plot_mat(u;name=string(oc.invε,"_",I,"_te_",round(Es0[I],digits=3),"_kin_",round(kin,digits=3)),dir="plots_wavefuns")
            end
        end
    end
end

function test_sym_R_sols_ex(ψs,oc,St_sol)
    @assert oc.e == "ex"
    @assert norm(oc.kred) < 1e-10 && abs(oc.V_intensity) > 1e-10 # to avoid loop
    be = get_basis_exact(oc,St_sol)
    p_exact = Params(2)
    p_exact.type = "mine_2d"
    p_exact.basis = be
    p_exact.us_K = ψs
    i = get_i_dmode(oc,St_sol)
    u1, u2 = get_u1_u2_in_good_basis(ψs[i], ψs[i+1], p_exact)
    ψs[i], ψs[i+1] = u1, u2
    test_rot_sym_u(u1,cis(2π/3),be;name="ψ_ex_1")
    test_rot_sym_u(u2,cis(4π/3),be;name="ψ_ex_2")
end

function compute_sol_ex(oc,St_sol::Stock)
    p = get_micro_p(oc,St_sol)
    be = get_basis_exact(oc,St_sol)
    ##
    Es, ψs = solve_exact(oc,St_sol;k_exact=true)#!only_Es)
    check_at_least_two_positive(Es)
    @assert norm(imag.(Es)) < 1e-10
    Es = real.(Es)
    if MAKE_TESTS && norm(oc.kred) < 1e-10 && abs(oc.V_intensity) > 1e-10
        test_sym_R_sols_ex(ψs,oc,St_sol) # can create loop
    end
    # i_center = get_i_dmode(oc,St_sol)
    # plot_Es_around(Es,i_center; amplitude=5)
    # proportions_fourier_modes(ψ,2,ifo.basis_exact;err=1e-4,name="exact_"*string(round(ε,digits=3)))
    ψs_reshaped = [myreshape(ψ_lin, be) for ψ_lin in ψs]
    Es, ψs_reshaped
end
