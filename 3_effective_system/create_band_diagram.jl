include("../1_common_graphene_manips/common_graphene.jl")

include("../lib/band_diagram.jl")
include("../lib/pt/asymptotic_behavior.jl")
include("../lib/rbm/test_rbm_pt.jl")

include("../2_monolayer_derivatives/graphene_derivatives.jl")

include("../3_effective_system/stock_module_solutions.jl")
include("../3_effective_system/compute_i_dmode.jl")
include("../3_effective_system/solve_effective.jl")
include("../3_effective_system/solve_exact.jl")
include("../3_effective_system/misc.jl")

include("../lib/stocks/stock_module.jl")
#################### Second step : compute the effective potentials 𝕍, 𝕎, 𝔸, etc. Very rapid

mean(l) = sum(l)/length(l)
path_pics_root = "/home/louis/Documents/tentatives/1-improved_graphene/pics/" # choose a convenient path for the pictures

# Computes effective potentials, makes plots, gives some information
function compute_bands(oc,St_sol;res_bands=10)
    px("Starts computing ",oc.e," band")
    send = true
    zoomed = false

    # p, ifo = generate_parameters(V_intensity=oc.V_intensity, order=oc.order, n_zeroth=oc.n_zeroth, k_dep=oc.k_dep,schur=oc.schur,V_kind=oc.V_kind,precision=oc.precision)
    bd = BandDiagram()
    bd.resolution_bands = res_bands
    # ifo.no_V = false # <---------------------------
    
    Γ, κ, K2, Γ2, M = build_Ks()
    zoom = 0.01
    K1p = zoom*κ; Mp = zoom*M
    # K1p = κ; Mp = M
    Klist = [κ,Γ,M]
    # Klist = [-κ,Γ,κ,2*κ]
    if zoomed
        Klist = [K1p,Γ,Mp]
    end
    Knames = ["κ","Γ","M"]
    # Knames = ["-κ","Γ","κ","2K_1"]
    init_klist_knames(Klist,Knames,bd)
    p = get_micro_p(oc,St_sol)
    add_lattice_vectors(p.basis.a1, p.basis.a2, bd)
    bd.k2lin = kred -> k2lin_ef(kred,oc,St_sol)

    # Choose nbands
    if oc.e == "ef"
        Nf = get_Nfull(oc,St_sol)
        bd.n_bands = Nf
    else
        bd.n_bands = get_N_basis_exact(oc,St_sol) # in stock_module_solution
        bd.n_bands = get_Nfull_exact(oc,St_sol)
        bd.n_bands = N_eigenmodes_exact(oc.invε,p)
    end
    function Hop(kred)
        new_oc = set_kred_to_oc(kred,oc)
        Es = get_Es(new_oc,St_sol)
        Es, ψs = X("sol",oc_to_list(new_oc),St_sol)
        # if oc.e == "ex" # plot exact 
            # plot_ψs_ex(new_oc.kred,Es,ψs,new_oc,St_sol)
            # be = get_basis_exact(oc,St_sol)
            # for u in ψs
                # test_nΩ_periodicity(u,oc.invε,be)
            # end
        # end
        #
        # oc_k0 = change_oc([0,0],"kred",oc)
        # Erf = get_Es(oc_k0,St_sol)[1]
        # Es .-= Erf
        Es
    end
    σ = spectrum_on_a_path(Hop,bd;print_progress=true,method="k_to_E")
    px("Ended computing ",oc.e," band")
    σ, bd
end

function plot_diag(oc,St_sol; only_one=false,more_name="",energy_amplitude=0.5, path = "",res_bands=10, print_N_macro=false, plot_yticks=false,finalized_plot=false)
    @assert only_one ∈ [false,"ex","ef"]
    print_oc(oc)
    σ_ef, σ_ex, bd, oc_ex, oc_ef = "", "", "", "", ""
    if only_one=="ex" || only_one==false 
        # σ_ex, bd, ifo_ex = X("σ_ex",dia_to_list(oc;exact=true),St_sol) # needs to be computed first because there is additional info in ifo and bd of effective
        
        oc_ex = set_e_to_oc("ex",oc)
        σ_ex, bd = compute_bands(oc_ex,St_sol;res_bands=res_bands)
        # dia = list_to_dia(coords)
        # compute_bands(dia,"effective")
    end
    if only_one=="ef" || only_one==false 
    # σ_ef, bd, ifo = X("σ_ef",dia_to_list(oc),St_sol)
    # σ_ex, bd, ifo_ex = compute_bands(oc,St_sol)
        oc_ef = set_e_to_oc("ef",oc)
        σ_ef, bd = compute_bands(oc_ef,St_sol;res_bands=res_bands)
    end
    # Ns = get_N_states(oc_ef,St_sol;total=false)
    bd.energy_center = 0
    bd.energy_amplitude = energy_amplitude
    bd.plot_all_bands = false
    bd.indicate_number_of_bands = true
    bd.graph_names = only_one==false ? ["ef","ex"] : [only_one]
    bd.plot_yticks = plot_yticks 
    # indicates i_mode
    i_modes_s = "i_modes "
    G = only_one==false ? 2 : 1
    for g=1:G
        e = bd.graph_names[g]
        noc = set_e_to_oc(e,oc)
        i_mode = get_i_dmode(noc,St_sol)
        i_modes_s = i_modes_s*e*"="*string(i_mode)*" "
    end
    # builds displayed info
    bd.displayed_info = [string("invε=",oc.invε,", ℓ=",oc.order) , string("N_macro=",oc.N_macro) , "V="*string(oc.V_intensity," kind=",oc.V_kind) , oc.k_dep , i_modes_s ]
    if oc.schur
        push!(bd.displayed_info,"Schur")
    end
    if finalized_plot
        bd.displayed_info = []
        bd.numbering_bands = false
    else
        bd.numbering_bands = true
    end
    resolution=(700,1000)
    # Save
    str_only_eff = only_one==false ? "" : "_"*only_one
    # str_no_V = ifo.no_V ? "_noV" : ""
    str_schur = oc.schur ? "_schur" : ""
    str_N_macro = print_N_macro ? "_Nmacro_"*string(oc.N_macro) : ""
    name = string("invep_",oc.invε,"_order_",oc.order,"_n0_",oc.n_zeroth,"_V_",oc.V_intensity,"_kind_",oc.V_kind,"_",oc.k_dep,str_only_eff,str_schur,str_N_macro,more_name)
    linestyles = [:dash,:dot]
    # px("Starts plotting bands")
    bd.root_path = "improved_graphene/plots_bands/"*string(oc.invε)*"/" #"./plots_bands/"*string(oc.invε)*"/"
    if path!=""
        bd.full_path = path
    end
    # if true
        # bd.root_path = [bd.root_path,"../../../tentatives/1-improved_graphene/pics/"]
    # end
    # if oc.out_name != nothing
        # name = oc.out_name
    # end
    bands = only_one==false ? [σ_ef, σ_ex] : (only_one=="ef" ? [σ_ef] : [σ_ex])
    colors = ["black","blue"]
    if only_one=="ef" colors = ["black"] end
    if only_one=="ex" colors = ["blue"] end
    plot_band_diagram(bands, "all", bd; post_name=name, extension="pdf", colors=colors, linewidth=0.25, linestyles=linestyles, resolution=resolution)
end

# creates the band diagram at ε=0 and V=0
function neutral_band_ε0_v0()
    St = create_stock_sol()
    oc = OneComp()
    # oc.N_macro = 5
    oc.precision = "precise"
    # oc.k_dep = "k_independent"
    oc.V_kind = "g"
    oc.invε = 1
    oc.V_intensity = 0
    # oc.order = 2
    res_bands = 40
    energy_amplitude = 1.6
    # oc.e = "ef"
    # i = get_i_dmode(oc,St)
    # px("i_dmode is ",i)
    plot_diag(oc,St;energy_amplitude=energy_amplitude,res_bands=res_bands,finalized_plot=true,only_one="ex",plot_yticks=true)
end

# Converge in N_macro : (ε,Nm) : (4,11(rather 13 for one upper))
function one_comp()
    St = create_stock_sol()
    oc = OneComp()
    # only_one = "ex"
    # only_one = "ef"
    # oc.N_macro = 10
    oc.precision = "precise"
    # oc.k_dep = "k_independent"
    oc.V_kind = "g"
    oc.invε = 1
    oc.V_intensity = 0
    # oc.order = 2
    res_bands = 20
    energy_amplitude = 0.79
    oc.e = "ex"
    # i = get_i_dmode(oc,St)
    # px("i_dmode is ",i)
    plot_diag(oc,St;energy_amplitude=energy_amplitude,res_bands=res_bands,finalized_plot=true)
end

function make_study()
    St = create_stock_sol()
    oc = OneComp()
    # only_one = "ex"
    # only_one = "ef"
    # only_one = false
    oc.N_macro = 5
    # oc.precision = "precise"
    invεs = [8,11,14,17,20]
    invεs = [11,19]
    invεs = [1,2,3,4,5]
    invεs = collect(6:15)
    invεs = [11,14,20]
    Vs = [0]
    schurs = [false]
    orders = [0]
    oc.V_kind = "g"
    energy_amplitude = 0.79
    res_bands = 5
    for invε in invεs, V in Vs, order in orders, schur in schurs
        # oc.V_kind = "rg"
        # oc.k_dep = "k_independent"
        for kdep in kdo(order)
            oc.invε = invε
            oc.V_intensity = V
            oc.schur = schur
            oc.order = order
            oc.k_dep = kdep
            plot_diag(oc,St;energy_amplitude=energy_amplitude,res_bands=res_bands,finalized_plot=false)
        end
    end
end

function study_ε(St0=0;short=true)
    St = St0==0 ? create_stock_sol() : St0
    oc = OneComp()
    # only_one = "ex"
    # only_one = "ef"
    only_one = false
    oc.precision = "precise"
    # oc.precision = "normal"
    invεs = [1,4,7,10,13,16]
    oc.order = 0
    oc.V_kind = "g"
    oc.k_dep = "k_independent"
    oc.V_intensity = 6
    # Plot options
    energy_amplitude = 0.79
    path = path_pics_root*"study_ep/"
    res_bands = 7
    N_macros = [2]
    if short
        oc.precision = "normal"
        # invεs = [1,5,7,11]
        res_bands = 5
    end
    for invε in invεs, N_macro in N_macros
        oc.N_macro = N_macro
        oc.invε = invε
        plot_diag(oc,St;print_N_macro=true,energy_amplitude=energy_amplitude,path=path,res_bands=res_bands,plot_yticks=invε==invεs[1],finalized_plot=true,only_one=only_one)
    end
    St
end

function study_ε_fixed_ν(St0=0;short=true)
    St = St0==0 ? create_stock_sol() : St0
    oc = OneComp()
    # only_one = "ex"
    only_one = "ef"
    oc.precision = "precise"
    # oc.precision = "normal"
    invεs = [7,10,19,31,40,49]
    oc.order = 1
    oc.V_kind = "g"
    oc.k_dep = "k_independent"
    oc.V_intensity = 0
    # Plot options
    energy_amplitude = 0.79
    path = path_pics_root*"study_ep_fixed_nu/"
    res_bands = 30
    N_macros = [21]
    if short
        oc.precision = "normal"
        res_bands = 4
    end
    for invε in invεs, N_macro in N_macros
        oc.invε = invε
        oc.N_macro = N_macro
        plot_diag(oc,St;print_N_macro=true,energy_amplitude=energy_amplitude,path=path,res_bands=res_bands,plot_yticks=invε==invεs[1],finalized_plot=true,only_one=only_one)
    end
    St
end

function study_ℓ(St0=0;short=true)
    St = St0==0 ? create_stock_sol() : St0
    oc = OneComp()
    oc.N_macro = 5
    oc.precision = "precise"
    oc.invε = 7
    orders = [0,1,2]
    vs = [6]
    oc.V_kind = "ng"
    # Plot options
    energy_amplitude = 0.79
    path = path_pics_root*"study_ell_kdep/"
    res_bands = 10
    if short
        oc.invε = 2
        res_bands = 7
        oc.precision = "normal"
    end
    for v in vs
        oc.V_intensity = v
        for order in orders
            for kdep in kdo(order)
                oc.k_dep = kdep
                oc.order = order
                oc.n_zeroth = 2
                plot_diag(oc,St;energy_amplitude=energy_amplitude,path=path,res_bands=res_bands,plot_yticks=order==0,finalized_plot=true)
            end
        end
        oc.n_zeroth = 6 # <--- excited states
        oc.order=0; oc.k_dep = "k_independent"; 
        plot_diag(oc,St;energy_amplitude=energy_amplitude,path=path,res_bands=res_bands,plot_yticks=false,finalized_plot=true)
    end
    St
end

function study_v(St0=0;short=true)
    St = St0==0 ? create_stock_sol() : St0
    oc = OneComp()
    oc.N_macro = 5
    oc.precision = "precise"
    oc.invε = 7
    vs = [0,1,5,15,30]
    vs = [0,1,3,4,5,6,7,8,9,10]
    oc.k_dep = "k_independent"
    # Plot options
    energy_amplitude = 0.79
    path = path_pics_root*"study_v/"
    res_bands = 30
    if short
        oc.invε = 2
        oc.precision = "normal"
        res_bands = 5
    end
    for v in vs
        oc.V_intensity = v
        oc.order = 2; oc.n_zeroth=2
        plot_diag(oc,St;energy_amplitude=energy_amplitude,path=path,res_bands=res_bands,plot_yticks=v==0,finalized_plot=true)
        # oc.order=0; oc.n_zeroth=6
        # plot_diag(oc,St;energy_amplitude=energy_amplitude,path=path,res_bands=res_bands,plot_yticks=j==0,finalized_plot=true)
    end
    St
end

function study_N(St0=0;short=true)
    St = St0==0 ? create_stock_sol() : St0
    oc = OneComp()
    N_macros = [3,5,7,9,13,17]
    oc.k_dep = "k_independent"
    oc.V_intensity = 0
    oc.invε = 7
    oc.precision = "precise"
    # Plot options
    energy_amplitude = 0.79
    path = path_pics_root*"study_N_macro/"
    res_bands = 50
    orders = [0,1,2]
    if short
        N_macros = [3,5,7]
        orders = [0]
        oc.invε = 2
        res_bands = 5
        oc.precision = "normal"
    end
    for N in N_macros, order in orders
        oc.order = order
        oc.N_macro = N
        plot_diag(oc,St;energy_amplitude=energy_amplitude,path=path,res_bands=res_bands,print_N_macro=true,plot_yticks=N==N_macros[1],finalized_plot=true)
    end
    St
end

short = false
# make_study()
St = 0
# St = study_ε(St,short=short)
# study_ε_fixed_ν(St;short=short)
# St = study_v(St,short=short)
# St = study_ℓ(St,short=short)
# St = study_N(St,short=short)
# St
neutral_band_ε0_v0()
# one_comp()
nothing
