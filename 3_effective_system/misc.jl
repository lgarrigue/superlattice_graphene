function get_MK()
    K = -[1/3,1/3]
    K1 = [-2,1]/3; K2 = [-1,2]/3
    Γ = [0.0,0.0]
    Γ2 = 2*K1-K2
    M = Γ2/2
    M, K
end

function build_Ks()
    Γ = [0.0,0.0]
    K1 = [-2,1]/3; K2 = [-1,2]/3
    Γ2 = 2*K1-K2
    M = Γ2/2
    Γ, K1, K2, Γ2, M
end

function modify_oc(λ,var,oc; V_base=0)
    if var=="N_macro"
        oc.N_macro = λ
    elseif var=="kred"
        oc.kred = λ
    elseif var=="invep"
        oc.invε = λ
    elseif var=="lambda_k"
        oc.kred = λ*K_DIR
    elseif var=="V"
        oc.V_intensity = λ
    elseif var=="V_and_ep"
        oc.V_intensity = V_base/λ
        oc.invε = λ
    end
    oc
end

function get_ψs_around(n,oc,St_sol; amplitude = 5)
    i = get_i_dmode(oc,St_sol)
    n0 = max(1,i-amplitude)
    N = get_Nfull(oc,St_sol)
    if oc.e=="ex"
        # N = get_Nfull_exact(oc,St_sol)
        p_micro = get_micro_p(oc,St_sol)
        N = N_eigenmodes_exact(oc.invε,p_micro)
    end
    n1 = min(N,i+amplitude)
    # px("n0 ",n0," n1 ",n1," N ",N)
    l = collect(n0:n1)
    ψs = []
    for p=n0:n1
        ψ = get_ψ(p,oc,St_sol)
        push!(ψs,ψ)
    end
    ψs, l
end

