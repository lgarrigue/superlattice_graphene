include("../../lib/plot/plot_matrix.jl")
include("../../lib/spectral_stuff.jl")

mutable struct State
    dependent_on_k
    data
    function State()
        sta = new()
        sta
    end
end

function get_state(sta,q)
    sta.dependent_on_k ? sta.data(q) : sta.data
end

function all_states_are_k_indep(ms)
    for sta in ms.states
        if sta.dependent_on_k
            return false
        end
    end
    true
end

# Stastes of the microscopic level
mutable struct MicroStates
    states
    all_states_k_indep
    previous_k_direction
    de
    schur
    function MicroStates(de)
        ms = new()
        ms.states = []
        # to avoid singularity when k=0
        ms.previous_k_direction = [1,1] 
        ms.de = de
        ms.schur = false
        ms
    end
end

order(ms) = ms.de.ℓ

function add_state_indep_k(ξ,ms)
    sta = State()
    sta.dependent_on_k = false
    sta.data = ξ
    push!(ms.states,sta)
end

function add_states_indep_k(ξs,ms)
    for (k,v) in ξs
        add_state_indep_k(v,ms)
    end
end

# add a state which is constant in k
function add_all_states_indep_k(ms)
    l = all_states(ms.de)
    # if gram
        # l = gram_schmidt_family(l)
    # end
    for ξ in l
        add_state_indep_k(ξ,ms)
    end
end

function add_states_from_fun(f,ms)
    for a=1:2
        sta = State()
        sta.dependent_on_k = true
        function g(kcd)
            # k = kcd
            # ms.previous_k_direction = kcd
            # if norm(kcd) < 1e-14
                # k = ms.previous_k_direction
                # k = -[1,0]
            # end
            k = kcd
            if norm(kcd) > 1e-14
                ms.previous_k_direction = kcd
            else
                k = ms.previous_k_direction
                ms.previous_k_direction *= -1
                if ms.previous_k_direction == [1,1] 
                    k = 1e-5*ms.previous_k_direction
                end
            end
            f(a,k,ms.de)
        end
        sta.data = g
        push!(ms.states,sta)
    end
end

# add a state which depends on the direction of k
function add_all_states_dep_k(ms)
    add_states_indep_k(ms.de.ξ0,ms)
    ℓ = order(ms)
    if ℓ ≥ 1
        add_states_from_fun(first_order_in_direction_k,ms)
    end
    if ℓ ≥ 2
        add_states_from_fun(second_order_in_direction_k,ms)
    end
end

function init_ms(de,k_dep)
    exact_eigenmodes = false
    ms = MicroStates(de)
    b = k_dep=="k_independent"
    ms.all_states_k_indep = b
    if b
        add_all_states_indep_k(ms)
    else
        add_all_states_dep_k(ms)
    end
    ms
end

######### Add R ∇v ϕ

function R_v_ϕ(p;n_scale=1,v_dir=p.v_dir)
    b = p.basis
    ξ0 = zeroth_order(p)
    v_fc = myfft(v_dir,b.N)
    if n_scale != 1
        v_fc = scale_2d(v_fc, n_scale, b, b) - v_fc
    end
    v_ϕ = [mul(myreshape(v_fc,b),ξ0[a],b) for a=1:2]
    [apply_pseudo_inverse(v_ϕ[a],p) for a=1:2]
end

function R_∇v_ϕ(p)
    b = p.basis
    ξ0 = zeroth_order(p)
    v_fc = myfft(p.v_dir,b.N)
    ∇v = ∇_k_mine(v_fc,b;kred=[0,0])
    ∇v_ϕ = [mul(myreshape(∇v[i],b),ξ0[a],b) for i=1:2, a=1:2]
    [apply_pseudo_inverse(∇v_ϕ[i,a],p) for i=1:2, a=1:2]
end

function R_dv_ϕ(p)
    v_fc = myfft(p.v_dir,p.basis.N)
    ∇v = ∇_k_mine(v_fc,p.basis;kred=[0,0])
    [apply_pseudo_inverse(∇v[xy],p) for xy=1:2]
end

function add_der_v(ms,p;nabla=false,v_dir=p.v_dir)
    f = nabla ? R_∇v_ϕ(p) : R_v_ϕ(p;v_dir=v_dir)
    for ξ in f
        add_state_indep_k(ξ,ms)
    end
end

# Gives the proportions of each state in one global state
function proportions_states(ψ,oc,St_sol)
    Ns = get_N_states(oc,St_sol)
    Nf = get_Nfull(oc,St_sol)
    props = zeros(Ns)
    for j=1:Nf
        (n,k) = i2m_and_k(j,oc,St_sol)
        props[n] += abs(ψ[j])
    end
    props
end

# Gives the proportions of each Fourier mode
function proportions_fourier_modes_α(αms,oc,St_sol)
    ψ = 0*copy(αms[1])
    for m=1:length(αms)
        ψ .+= abs.(αms[m])
    end
    proportions_fourier_modes(ψ,2,get_basis_macro(oc,St_sol))
end

# creates the preconditionner for LOBPCG
function k2lin_ef(kred_base,oc,St_sol) 
    l = []
    Nf = get_Nfull(oc,St_sol)
    for i=1:Nf
        (m,k_red) = i2m_and_k(i,oc,St_sol)
        k_cart = kred2kcart_2d(k_red,get_basis_macro(oc,St_sol);shift=kred_base)
        k2 = norm(k_cart)^2
        push!(l,k2)
    end
    l
end

i2m_and_k(i,oc,St_sol) = get_map(oc,St_sol)[i]

# Direction of kcart
# d for direction, r for rounded
# we round to ensure that it is the same that is taken in dict entries, this is not used for exact purposes
function kcart_dr(kred,oc,St_sol) 
    ms = get_ms(oc,St_sol)
    # px("COMPUTES DRI")
    kcart = kred2kcart_2d(kred,get_basis_macro(oc,St_sol))
    if norm(kcart)<1e-13 || ms.all_states_k_indep
        return [0,0] # was at [1,1] before
    end
    x = kcart/norm(kcart) 
    x = round.(x,digits=12)
    return x
end

function heart_sing_op(m1,m2,kred,kred_global,oc,St_sol)
    St_ms = get_St_ms_coefs(oc,St_sol)
    kcd = kcart_dr(kred_global,oc,St_sol)
    x = X("ξ(h-E)ξ", (m1,m2,kcd), St_ms)
    # px("ξ(h-E)ξ",abs(x))
    # @assert abs(x)<1e-13
    x
end

function heart_reg_op(m1,m2,kred,kred_global,oc,St_sol)
    St_ms = get_St_ms_coefs(oc,St_sol)
    kcd = kcart_dr(kred_global,oc,St_sol)
    a,b = X("ξ(-i∇)ξ", (m1,m2,kcd), St_ms)
    kr = kred + kred_global
    k = kred2kcart_2d(kr,get_basis_macro(oc,St_sol))
    # @assert (abs(a)<1e-13 || abs(a-0.38766)<1e-3)
    k[1]*a + k[2]*b
end

function heart_last_op(m1,m2,kred,kred_global,oc,St_sol)
    St_ms = get_St_ms_coefs(oc,St_sol)
    kcd = kcart_dr(kred_global,oc,St_sol)
    ps = X("ξξ", (m1,m2,kcd), St_ms)
    kr = kred + kred_global
    k = kred2kcart_2d(kr,get_basis_macro(oc,St_sol))
    0.5*norm(k)^2*ps
end

function mass_op(m1,m2,kred,kred_global,oc,St_sol)
    St_ms = get_St_ms_coefs(oc,St_sol)
    kcd = kcart_dr(kred_global,oc,St_sol)
    x = X("ξξ", (m1,m2,kcd), St_ms)
    # px("x= ",x," m1 m2 kcd ",m1," ",m2," ",kcd)
    x
end

function heart_pot_op(m1,m2,kreds,kred_global,oc,St_sol)
    macro_b = get_basis_macro(oc,St_sol)
    St_ms = get_St_ms_coefs(oc,St_sol)
    k1, k2 = kreds
    kcd = kcart_dr(kred_global,oc,St_sol)
    ps = X("ξξ", (m1,m2,kcd), St_ms)
    kr = k1-k2
    bm = get_basis_macro(oc,St_sol)
    ik = kred2ik_2d(kr, bm)
    # c = length(ifo.V_fc)/sqrt(bm.Vol)
    # l = length(ifo.V_fc)
    s = sqrt(bm.Vol)
    # px("l, s, c = ",l," ",s," ",c)
    if ik!=-1
        x = get_V_fc(oc,St_sol)[ik]*ps/s
        # if abs(x)>1e-14
            # px("x is ",x," ik ",ik," kred ",kr)
        # end
        return x
    else
        return 0
    end
end

function schur_op(m1,m2,kred,kred_global,oc,St_sol)
    St_ms = get_St_ms_coefs(oc,St_sol)
    kcd = kcart_dr(kred_global,oc,St_sol)
    kr = kred + kred_global
    k = kred2kcart_2d(kr,get_basis_macro(oc,St_sol))
    sum(-k[a]*k[b]*X("TMT",(a,b,kcd),St_ms)[m1,m2] for a=1:2, b=1:2)
end

# fourier_diag tels is it's diagonal action in Fourier space
function create_op(fun,oc,St_sol;fourier_diag=true)
    k_red_global = oc.kred
    Nf = get_Nfull(oc,St_sol)
    A = 0*rand(ComplexF64,Nf,Nf)
    for i=1:Nf
        # px("i = ",i)
        (m1,k1) = i2m_and_k(i,oc,St_sol)
        for j=1:Nf
            (m2,k2) = i2m_and_k(j,oc,St_sol)
            if fourier_diag
                if k1==k2
                    A[i,j] = fun(m1,m2, k1, k_red_global,oc,St_sol)
                end
            else
                A[i,j] = fun(m1,m2, (k1,k2), k_red_global,oc,St_sol)
            end
        end
    end
    A
end

function analyze_sing_op(kred,oc,St_sol)
    sing_op = create_op(kred,heart_sing_op,oc,St_sol)
end

# qred is the momentum shift
function ℍ(oc,St_sol;V_equal_0=false)
    ε = 1/oc.invε
    sing_op = create_op(heart_sing_op,oc,St_sol)
    reg_op  = create_op(heart_reg_op,oc, St_sol)
    last_op = create_op(heart_last_op,oc,St_sol)
    H = (1/ε)*sing_op + reg_op + ε*last_op
    
    if abs(oc.V_intensity) > 1e-14
        V_op = create_op(heart_pot_op,oc,St_sol; fourier_diag=false)
        H += V_op
    end
    if oc.schur
        so = create_op(schur_op,oc,St_sol)
        assert_herm(so,1e-10)
        px("REGLER LE BON SCLAING SUR SCHUR")
        H += ε*so
    end
    if MAKE_TESTS
        assert_herm(sing_op,1e-10)
        assert_herm(reg_op,1e-10)
        assert_herm(last_op,1e-10)
    end
    H
end

𝕊(oc,St_sol) = create_op(mass_op,oc,St_sol)

function print_macro_vec(ψ,oc,St_sol)
    Nf = get_Nfull(oc,St_sol)
    for i=1:Nf
        c = abs(ψ[i])
        if c > 0.00001
            (m,k) = i2m_and_k(i,oc,St_sol)
            px("m=",m," k=",k," |ψ|=",c)
        end
    end
end

conv(u,v) = FFTW.fft(FFTW.ifft(u) .* FFTW.ifft(v))

### Extracts the macroscopic eigenmodes
# From β (α linear) to α(εx)
function β_to_α(β,oc,St_sol)
    Nf = get_Nfull(oc,St_sol)
    macro_b = get_basis_macro(oc,St_sol)
    Ns = get_N_states(oc,St_sol)
    # Build α_{m}'s
    α = [0*copy(new_fun_2d(macro_b)) for m=1:Ns]
    for i=1:Nf
        (m,kred_α) = i2m_and_k(i,oc,St_sol)
        ik = kred2ik_2d(kred_α,macro_b)
        # px("ik=",ik," kred_α=",kred_α," Nf=",Nf)
        if ik!=-1
            α[m][ik] = β[i]
        else
            px("PROBLEM IN RECONSTRUCTING α") # maybe it's not a problem...
        end
    end
    α
end

function α_to_αψ(α,oc,St_sol)
    p = get_micro_p(oc,St_sol)
    mb = get_micro_basis(oc,St_sol)
    macro_b = get_basis_macro(oc,St_sol)
    be = get_basis_exact(oc,St_sol)
    St_ms = get_St_ms_coefs(oc,St_sol)
    Ns = get_N_states(oc,St_sol)
    ### Multiply with base functions. αψs are in fc
    αψ = 0*zeros(ComplexF64, be.dims)
    for m=1:Ns
        kcd = kcart_dr(oc.kred,oc,St_sol)
        ψ = X("ξ",(m,kcd),St_ms)
        ψ_fc = u_lin2fc(ψ,p)
        α_ = change_N_matrix(α[m],macro_b, be) 
        # if p.basis.dim==3
            # for ikz=1:p.basis.Nz
                # φ = ψ_fc[:,:,ikz]
                # ψ_scaled = scale_2d(ψ_fc, oc.invε, mb, be)
                # αψs[j][:,:,ikz] += conv(α,ψ_scaled) # the good one when things work
            # end
        # else
            # n_scale = ifo.no_V ? 1 : n
        ψ_scaled = scale_2d(ψ_fc, oc.invε, mb, be)
        αψ += conv(α_,ψ_scaled)
        # end
    end
    # Normalize
    no = norm(αψ) # is sometimes 0
    if no > 1e-10
        αψ /= no
    else
        if abs(oc.V_intensity) > 1e-14 # if V!=0, then αs should be non zero
            px("PROBLEM NORM in eff_pot, norm of αψ is ",no)
        end
    end
    αψ
end
   
function compute_sol_ef(oc,St_sol)
    t0 = time()
    H = ℍ(oc,St_sol)
    S = 𝕊(oc,St_sol)
    # px("Built matrices effective in ",time()-t0," sec")
    ### Solves the macroscopic eigenmodes
    tol = 1e-9
    if MAKE_TESTS px("IS SLOWED BY CHECK TESTS") end
    Es, βs = gen_eigenmodes(H,S; check_tests=MAKE_TESTS,tol=tol) # in spectral_stuff.jl, exact one
    # px("Solved effective after ",time()-t0," sec")
    # bm = get_basis_macro(oc,St_sol)
    # px("Infos effective : Ns=",get_N_states(oc,St_sol),", Nfull=",get_Nfull(oc,St_sol),", N_macro=",bm.N)
    # αψs, αs = build_effective_eigenvects(βs,oc,St_sol)
    check_at_least_two_positive(Es)
    @assert norm(imag.(Es)) < 1e-10
    Es = real.(Es)
    Es, βs #αψs
end
