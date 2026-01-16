# using CairoMakie
# include("../../lib/stock_module.jl")
# px = println
include("../3_effective_system/stock_module_coefs.jl")

# HARDCODED MACROS
a_cell = 5
i_state_micro = 1
MAKE_TESTS = false
COEF_v_micro = 10
#

function prec2N(precision;micro_basis=false,macro_basis=false)
    N = 9
    if precision=="no"
        N = 5
    elseif precision=="normal"
        N = 7
    elseif precision=="ultra"
        N = 11
    end
    if micro_basis
        nothing
    end
    # if macro_basis
        # return 9
    # end
    N
end

kdo(ℓ) = ℓ==0 ? ["k_independent"] : ["k_dependent","k_independent"]

# Contains data of one computation
mutable struct OneComp
    e
    invε
    V_intensity
    V_kind
    kred

    N_macro
    order
    n_zeroth
    k_dep
    schur
    precision

    function OneComp()
        oc = new()
        oc.e = "ef"
        oc.invε=1
        oc.V_intensity=0
        oc.V_kind="g"
        oc.kred=[0,0]

        oc.N_macro=7
        oc.order=0
        oc.n_zeroth=2
        oc.k_dep="k_independent"
        oc.schur=false
        oc.precision = "normal"
        oc
    end
end

function print_oc(oc)
    px("OC ",oc.e," 1/ε=",oc.invε," precision ",oc.precision," N_macro=",oc.N_macro," V=",oc.V_intensity," kind=",oc.V_kind," k=",oc.kred," | ℓ=",oc.order," n_zeroth=",oc.n_zeroth," ",oc.k_dep," schur? ",oc.schur)
end

function oc_to_list(oc)
    l = [oc.e,oc.invε,oc.V_intensity,oc.V_kind,oc.kred,oc.precision]
    if oc.e=="ex"
        return l
    end
    vcat(l,[oc.N_macro,oc.order,oc.n_zeroth,oc.k_dep,oc.schur])
end

function list_to_oc(l)
    oc = OneComp()
    # Exact and effective system parameters
    oc.e = l[1]
    oc.invε = l[2]
    oc.V_intensity = l[3]
    oc.V_kind = l[4]
    oc.kred = l[5]
    oc.precision = l[6]
    if oc.e == "ex"
        return oc
    end
    # Only effective system parameters
    oc.N_macro = l[7]
    oc.order = l[8]
    oc.n_zeroth = l[9]
    oc.k_dep = l[10]
    oc.schur = l[11]
    oc
end

variables_link = Dict()
variables_link["kred"] = 5
variables_link["k_dep"] = 10
variables_link["invε"] = 2
variables_link["V"] = 3
variables_link["e"] = 1
variables_link["n_zeroth"] = 9
variables_link["order"] = 8

function change_oc(λ,var,oc)
    noc_l = oc_to_list(oc)
    noc_l[variables_link[var]] = λ
    noc = list_to_oc(noc_l)
    noc
end

function set_kred_to_oc(kred,oc)
    noc = oc; noc.kred = kred
    noc
end

function set_e_to_oc(e,oc)
    noc = oc; noc.e = e
    noc
end


function create_stock_sol()
    pS = Dict()
    # pS["ifo"] = ifo

    function compute_sol(coords,St_sol::Stock)
        oc = list_to_oc(coords)
        if oc.e=="ex"
            Es, ψs_reshaped = compute_sol_ex(oc,St_sol)
            return Es, ψs_reshaped
        end
        compute_sol_ef(oc,St_sol)
    end
    function compute_p_micro(coords,St_sol) # p
        precision = coords
        ################# Build micro
        N_micro = prec2N(precision) # <-------------- increase if precision is not enough. Take 5 if really want to make it quick, but does not really improve
        # N_micro = 5 # <-------------- increase if precision is not enough, should be 20 ?
        v_intensity = COEF_v_micro # positive or negative doesn't matter, we are in Fourier space with oscillating functions !
        p = generate_some_2d_param(N_micro, v_intensity; i_state=i_state_micro, a=a_cell)
        # px("NNNNNNNNNNNNNN ",mean(abs.(p.v_dir)))
        p
    end
    function compute_basis_macro(coords,St_sol)
        # @assert length(coords) = 1
        N_macro = coords
        basis_macro = graphene_basis_2d(N_macro; a=a_cell)
        basis_macro
    end
    function compute_basis_exact(coords,St_sol::Stock)
        invε, precision = coords
        N_micro_basis = prec2N(precision)
        N = N_basis_exact(invε,N_micro_basis;precision=precision)
        basis_exact = graphene_basis_2d(N;a=a_cell)
        basis_exact 
    end
    function compute_V_dir(coords,St_sol)
        V_intensity, V_kind, N_macro = coords
        basis_macro = X("basis_macro",N_macro,St_sol)
        V_dir = ""
        R_four_2d = matrix_rot_red(-2π/3,basis_macro)
        if V_kind=="g" # graphene
            V_dir = generate_some_sym_potential_dir(V_intensity, basis_macro, first_k=[1,0])
        elseif V_kind=="rg" # rotated graphene
            fk = [1,1]
            V_dir = generate_some_sym_potential_dir(V_intensity, basis_macro, first_k=fk)
            # kreds = [(R_four_2d^i)*fk for i=0:3]
            # px("KREDS ",kreds)
        elseif V_kind=="ng" # not graphene
            # fk = [-1,1]
            # kreds = [(R_four_2d^i)*fk for i=0:3]
            kreds = [[1,0],[0,2],[0,1]]
            kreds = vcat(kreds,-kreds)
            # px("##################### NG POT\n\n")
            V_dir = generate_some_sym_potential_dir(V_intensity, basis_macro, mykreds=kreds)
        else
            px("PROBLEM IN V_kind")
            @assert false
        end
        V_dir
        # ifo.V_dir = V_dir
        # ifo.V_fc = myfft(V_dir,ifo.basis_macro.Vol)
        # ifo.V_intensity = V_intensity
        # V_fc = change_N_matrix(myfft(V_dir,ifo.basis_macro.N), ifo.basis_macro, p.basis)
        # stock = init_stock_micro_coefs(p; ℓ=order, typ=k_dep, V_dir=myifft(V_fc,p.basis.N), n_zeroth=n_zeroth, exact_eigenmodes=exact_eigenmodes, kred_exact=kred_exact) # yavait pas besoin de V_dir la dedans je pense
    end
    function compute_V_fc(coords,St_sol)
        V_intensity, V_kind, N_macro = coords
        basis_macro = X("basis_macro",N_macro,St_sol)
        V_dir = X("V_dir",coords,St_sol)
        myfft(V_dir,basis_macro.Vol)
    end
    function compute_ef_map(coords,St_sol::Stock)
        order, n_zeroth, k_dep, N_macro = coords
        basis_macro = X("basis_macro",N_macro,St_sol)
        St_ms = X("St_ms_coefs",coords,St_sol)
        ms = St_ms.p["ms"]
        Ns = length(ms.states)
        # px("Ns ",Ns," N ",basis_macro.N)
        # px("Number of micro states ",Ns)
        map = [(m, basis_macro.kred_grid_2d[i,j]) for m=1:Ns, i=1:basis_macro.N, j=1:basis_macro.N] # if basis_macro.kred_grid_2d is not defined, it's because basis_macro is not well defined yet
        map
    end
    function compute_St_ms_coefs(coords,St_sol::Stock)
        order, n_zeroth, k_dep, precision = coords
        p = X("p_micro",precision,St_sol)
        exact_eigenmodes = false
        kred_exact=[0,0]
        de = Derivatives(order,p; n_zeroth=n_zeroth, di_triplet=3,exact_eigenmodes=exact_eigenmodes,kred_exact=kred_exact)
        ms = init_ms(de,k_dep)
        gram = true
        if order==2 && k_dep=="k_independent"
            gram = true
        end
        create_stock_micro_coefs(ms, p; gram=gram)
    end

    kcf = [("sol",compute_sol),
           ("i_dmode",compute_i_dmode),
           ("p_micro",compute_p_micro),
           ("basis_macro",compute_basis_macro),
           ("basis_exact",compute_basis_exact),
           ("V_dir",compute_V_dir),
           ("V_fc",compute_V_fc),
           ("v_ex_dir",compute_exact_potential),
           ("ef_map",compute_ef_map),
           ("St_ms_coefs",compute_St_ms_coefs)]
    St_sol = Stock(kcf,pS)
    St_sol
end

function get_Es_unshifted(oc,St_sol) # Es which are not shifted at k=0 energy
    X("sol",oc_to_list(oc),St_sol)[1]
end

function get_Es(oc,St_sol)
    Es = get_Es_unshifted(oc,St_sol)
    i = get_i_dmode(oc,St_sol)
    oc_k0 = change_oc([0,0],"kred",oc)
    # oc_k0 = oc; oc_k0.kred = [0,0]
    Es_k0 = get_Es_unshifted(oc_k0,St_sol)
    Es .- Es_k0[i]
end

function get_ψ(n,oc,St_sol) # nth eigenstate
    γ = X("sol",oc_to_list(oc),St_sol)[2][n]
    if oc.e=="ex"
        return γ
    end
    α = β_to_α(γ,oc,St_sol)
    α_to_αψ(α,oc,St_sol)
end

function get_α_ef(n,oc,St_sol) # nth eigenstate
    β = X("sol",oc_to_list(oc),St_sol)[2][n]
    β_to_α(β,oc,St_sol)
end

function get_basis_macro(oc,St_sol)
    X("basis_macro",oc.N_macro,St_sol)
end

function get_micro_p(oc,St_sol)
    X("p_micro",oc.precision,St_sol)
end

function get_N_micro(oc,St_sol)
    X("p_micro",oc.precision,St_sol).basis.N
end

function get_micro_basis(oc,St_sol)
    X("p_micro",oc.precision,St_sol).basis
end

function get_V_dir(oc,St_sol)
    coords = oc.V_intensity, oc.V_kind, oc.N_macro
    X("V_dir",coords,St_sol)
end

function get_V_fc(oc,St_sol)
    coords = oc.V_intensity, oc.V_kind, oc.N_macro
    X("V_fc",coords,St_sol)
end

########### Effective gets

function get_N_states(oc,St_sol;total=false)
    if oc.schur && !total return 2 end
    ms = get_ms(oc,St_sol)
    length(ms.states)
end

function get_Nfull(oc,St_sol)
    map = get_map(oc,St_sol)
    # px("map is ")
    # display(map)
    length(map)
end

function get_map(oc,St_sol)
    coords = oc.order, oc.n_zeroth, oc.k_dep, oc.N_macro
    X("ef_map",coords,St_sol)
end

function get_St_ms_coefs(oc,St_sol)
    coords = oc.order, oc.n_zeroth, oc.k_dep, oc.precision
    X("St_ms_coefs",coords,St_sol)
end

function get_ms(oc,St_sol)
    St = get_St_ms_coefs(oc,St_sol)
    St.p["ms"]
end


########### Exact gets

function get_basis_exact(oc,St_sol)
    X("basis_exact",(oc.invε,oc.precision),St_sol)
end

function get_N_basis_exact(oc,St_sol)
    p = get_micro_p(oc,St_sol)
    N = N_basis_exact(oc.invε,p.basis.N)
    N^2
end

function get_Nfull_exact(oc,St_sol)
    be = get_basis_exact(oc,St_sol)
    be.N^2
end
