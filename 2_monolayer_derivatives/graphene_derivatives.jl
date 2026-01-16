include("../lib/linalg/linear_independence.jl")
# Derivatives with respect to k

mutable struct Derivatives
    n_zeroth # number of states in the zeroth order
    i_doublet # position of the Dirac doublet
    di_triplet # position of the singlet from the triplet corresponding to the targeted Dirac Bloch eigenfunctions
    params
    normalize
    ℓ
    ξ0
    ξ1
    noR_ξ1
    ξ2
    ξ2_deg
    ξ2_dir_k_prepa
    P
    function Derivatives(ℓ,p;n_zeroth=2, di_triplet=3,exact_eigenmodes=false,kred_exact=[0,0])
        de = new()
        de.n_zeroth = n_zeroth
        de.params = p
        de.i_doublet = 1 
        de.di_triplet = "" # to show that it does nothing
        de.normalize = false
        de.ℓ = ℓ
        init_Derivatives(de;exact_eigenmodes=exact_eigenmodes,kred_exact=kred_exact)
        de
    end
end

function all_ξ2(de)
    l = []
    li = [de.ξ2, de.ξ2_deg]
    li = [de.ξ2_dir_k_prepa]
    for list in li
        for ξ in list
            push!(l,ξ)
        end
    end
    # px("NUMBER OF VECTORS GIVEN BY all_ξ2 ",length(l))
    l
end

function init_Derivatives(de;exact_eigenmodes=false,kred_exact=[0,0])
    compute_zeroth_order(de;exact_eigenmodes=exact_eigenmodes,kred_exact=kred_exact)
    compute_P(de)
    compute_first_order(de)
    compute_second_order(de)
end

function all_states(de)
    l = []
    ds = [de.ξ0, de.ξ1, all_ξ2(de)]
    for j=1:1+de.ℓ
        for (k,v) in ds[j]
            push!(l, v)
        end
    end
    # test_linear_independence(l;tol=1e-4)
    l
end

function all_directional_states(dk,de)
    b = de.params.basis
    kc = kred2kcart_dim(dk,b)
    l = []
    function add(v)
        push!(l,v)
    end
    # add(de.ξ0)
    for a=1:2
        add(de.ξ0[a])
        if de.ℓ ≥ 1
            M1 = first_order_in_direction_k(a,kc,de)
            add(M1)
            if de.ℓ ≥ 2
                M2 = second_order_in_direction_k(a,kc,de)
                add(M2)
            end
        end
    end
    # test_linear_independence(l;tol=1e-4)
    l
end

# a ∈ [1,2] corresponds to u1 and u2
# Derivative with respect to k
mi∇(ψ,p) = (-im).*diff_op_K(ψ,"∇_K",p)

function R_mi∇(ψ,p)
    mi∇u = mi∇(ψ,p)
    [apply_pseudo_inverse(mi∇u[xy],p) for xy=1:2]
    # [∇u[xy] for xy=1:2]
end

function R2_mi∇(ψ,p)
    Rmi∇u = R_mi∇(ψ,p)
    [apply_pseudo_inverse(Rmi∇u[xy],p) for xy=1:2]
    # [∇u[xy] for xy=1:2]
end

# P2(u,v) = u*u' + v*v'

# N_triplet ∈ [2,-2]
function compute_zeroth_order(de;exact_eigenmodes=false,kred_exact=[0,0]) # if exact, then we consider the exact Bloch eigenstate
    d = Dict()
    us = []
    if exact_eigenmodes
        p = de.params
        K = p.basis.K_red_2d
        k = kred_exact+K
        _, us = diag_monolayer_at_k_mine(k,p.v_dir,p.basis.Nfull,p.basis;tol=1e-5)#,maxiter=100,full_diag=true,coef_kin=1)
        # the regular call is in my_2d_structs
        # px("############### istate ",p.i_state," V ",norm(p.v_dir)," l ",p.basis.Nfull," DIST ",dist_quo(us[1],u_K(1,p)))
        # u_K(1) is rotated so it's normal it's not equal to us[1]
        # px("############### ",norm(us[1]*us[1]' + us[2]*us[2]' - u_K(1,p)*u_K(1,p)' - u_K(2,p)*u_K(2,p)')) was 0 so it's ok
    end
    for i=1:de.n_zeroth
        I = i
        # I = i==3 ? de.di_triplet+1 : i # if we want the third excited state, give the triplet number, which could be below the doublet
        uk = u_K(I,de.params)
        if exact_eigenmodes
            uk = us[I]
        end
        d[i] = uk
    end
    de.ξ0 = d
end

function normalize_dict(d)
    for (k,v) in d
        d[k] /= norm(d[k])
    end
end

# ξ1[a,α] = R (-i∂α) u_a
function compute_first_order(de) # 4 states
    de.ξ1 = Dict()
    de.noR_ξ1 = Dict()
    for a=1:2, α=1:2
        de.ξ1[(a,α)] = R_mi∇(de.ξ0[a], de.params)[α]
        de.noR_ξ1[(a,α)] = mi∇(de.ξ0[a], de.params)[α]
    end
    if de.normalize
        normalize_dict(de.ξ1)
        normalize_dict(de.noR_ξ1)
    end
end

function compute_P(de) # the projection onto the two Bloch eigenstates
    i = de.i_doublet
    u1 = u_K(i,  de.params)
    u2 = u_K(i+1,de.params)
    de.P = u1*u1' + u2*u2'
end

# ξ2[a,β,α] = R (-i∂β) R (-i∂α) u_a
# ξ2_deg[a,α] = R^2 (-i∂α) u_a
function compute_second_order(de) # 8 + 4 states
    de.ξ2 = Dict()
    de.ξ2_deg = Dict()
    de.ξ2_dir_k_prepa = Dict()
    for a=1:2, α=1:2, β=1:2
        # R H1 R H1 ψ
        de.ξ2[(a,β,α)] = R_mi∇(de.ξ1[(a,α)], de.params)[β]
        de.ξ2_dir_k_prepa[(a,β,α)] = de.ξ2[(a,β,α)] - apply_pseudo_inverse(R_mi∇(de.P*mi∇(de.ξ0[a], de.params)[α], de.params)[β], de.params)
    end
    # R^2 H1 P H1 ψ
    for a=1:2, α=1:2
        de.ξ2_deg[(a,α)] = apply_pseudo_inverse(de.ξ1[(a,α)] ,de.params)
    end
    if de.normalize
        normalize_dict(de.ξ2)
        normalize_dict(de.ξ2_deg)
    end
end

function print_props(φs)
    Ns = length(φs)
    norms = [norm(φ) for φ in φs]
    su = [sum(norms[1:2])]
    if Ns ≥ 3
        push!(su,sum(norms[3:min(6,Ns)]))
    end
    if Ns ≥ 7
        push!(su,sum(norms[7:Ns]))
    end
    px("Weights ",su)
    su
end

# ξ1 = first_order(p;nor=false)
# a ∈ [1,2]
# kc = kcart
# Dk(u,kc,de) = sum(kc[α]*de.ξ1[(a,α)] for α=1:2)
first_order_in_direction_k(a,kc,de) = sum(kc[α]*de.ξ1[(a,α)] for α=1:2)
# R_first_order_in_direction_k(a,kc,de) = sum(kc[α]*de.ξ1[(a,α)] for α=1:2)
function second_order_in_direction_k(a,kc,de)
    sum(kc[α]*kc[β]*de.ξ2_dir_k_prepa[(a,α,β)] for α=1:2, β=1:2)
end

function M_first_order(a,b,de;R=false)
    M = zeros(ComplexF64,2,2)
    for α=1:2, β=1:2
        f = R ? de.ξ1 : de.noR_ξ1
        M[α,β] = f[(a,α)]⋅de.ξ1[(b,β)]
    end
    M
end

# < R (-i∂_α) w_a , (R) (-i∂_β) R (-i∂_γ) w_b > 
function F_second_order(de;R=false)
    F = Dict()
    for α=1:2, β=1:2, γ=1:2, a=1:2, b=1:2
        ψ = de.ξ1[(b,γ)]
        f = R ? R_mi∇ : mi∇
        φ = f(ψ,de.params)[β]
        F[(a,b),(α,β,γ)] = de.ξ1[(a,α)]⋅φ
    end
    F
end

function F_vecs(a,b,de;R=false)
    F = F_second_order(de;R=R)
    coords = get_coords(;dims=3)
    Fvec = zeros(ComplexF64,8)
    for A=1:8
        Fvec[A] = F[(a,b),coords[A]]
    end
    Fvec
end

function small(x;tol=1e-13)
    n = norm(x)
    b = n < tol
    if !b
        px("PROBLEM IN SMALL ",n)
    end
    return n, b
end

function print_first_order(de)
    for R in [false,true]
        px("####### First order, R=",R)
        D = [1 im; im -1]
        σ2 = [0 -im;
              im 0]

        M11 = M_first_order(1,1,de;R=R)
        M12 = M_first_order(1,2,de;R=R)
        M21 = M_first_order(2,1,de;R=R)
        M22 = M_first_order(2,2,de;R=R)

        t = M11[1,1]
        s = M11[2,1]/im
        r = M12[1,1]
        labs = ["r","t","s"]
        vals = [r,   t,  s]
        disp_vals(vals,labs)
        px("(t+s)/2 +r=",(t+s)/2 +r)

        px("Check tI + s σ2 = M11")
        small(t*I + s*σ2 - M11)
        px("Check r D = M12")
        small(r*D - M12)
        px("Check M12* = M21")
        small(M12' - M21)
        px("Check conj(M22) = M11")
        small(conj.(M22) - M11)
    end
end

# asserts that the numbers in l are all equal
function assert_equality(l;tol=1e-10)
    ref = l[1]
    for i=2:length(l)
        small(ref-l[i])
    end
end

function disp_vals(vals,labs;tol=1e-7)
    for i=1:length(labs)
        @assert abs(imag(vals[i]))<tol
        # px(labs[i],"= ",vals[i])
    end
    vals = real.(vals)
    for i=1:length(labs)
        px(labs[i],"= ",vals[i]," Hartree, ",vals[i]*hartree_to_ev," eV")
    end
end

function ul(a)
    if a==1
        return 2
    elseif a==2
        return 1
    end
    px("PROBLEM IN UL")
    -9
end

function theory(a,b,α,β,γ,χγ) # theoreticall values of F^{ab}_{αβγ}
    if (a,b) in [(1,1),(2,2)]
        ε1 = (α==β && β==γ) ? 1 : -1
        ε2 = a==2 ? -1 : 1
        bo = (α+β+γ)%2==0
        χ = bo ? ε2*χγ["χ2"] : χγ["χ1"]
        return ε1*χ
    else
        γ1 = χγ["γ1"]; γ2 = χγ["γ2"]; γ3 = χγ["γ3"]
        return 0
    end
end

function print_second_order(de;tol=1e-6)
    for R in [false,true]
        px("####### Second order, R=",R)
        F = F_second_order(de;R=R)

        χ1 = F[(1,1),(1,1,1)]
        χ2 = -im*F[(1,1),(2,2,2)]
        γ1 = F[(1,2),(1,2,2)]
        γ2 = F[(1,2),(2,1,2)]
        γ3 = F[(1,2),(2,2,1)]
        labs = ["χ1","χ2","γ1","γ2","γ3","γ1+γ2+γ3"]
        vals = [χ1,   χ2,  γ1,  γ2,  γ3,  γ1+γ2+γ3]
        # px(labs[i],"= ",vals[i])
        # end
        # vals = real.(vals)
        χγ = Dict()
        for i=1:length(labs)
            v = vals[i]
            @assert abs(imag(v))<tol
            rv = real(v)
            px(labs[i],"= ",rv," Hartree, ",rv*hartree_to_ev," eV")
            χγ[labs[i]] = rv
        end

        function x_to_str(x)
            s = "?"; tol = 1e-10; count = 0
            for j=1:length(vals)
                x_lab = labs[j]
                v = vals[j]
                if abs(v-x)<tol
                    s = labs[j]
                    count += 1
                elseif abs(v+x)<tol
                    s = "-"*labs[j]
                    count += 1
                elseif abs(x+im*v)<tol
                    s = "-i"*labs[j]
                    count += 1
                elseif abs(x-im*v)<tol
                    s = "i"*labs[j]
                    count += 1
                end
            end
            if count >= 2
                println("COUNT = ",count)
                # @assert false
            end
            return s
        end

        for (a,b) in [(1,1),(2,2),(1,2),(2,1)]
            for α=1:2, β=1:2, γ=1:2
                ex = F[(a,b),(α,β,γ)]
                s = x_to_str(ex)
                px(a,b,"  ",α,β,γ," ",s)
                ap = theory(a,b,α,β,γ,χγ)
                # if abs(ap-ex)>tol
                # px("NOT GOOD, is ",s," should be ",x_to_str(ap))
                # px("exact=",ex," approx=",ap," chi",bo ? "2" : "1")
                # end
            end
        end

        # Test symmetries from frakR and CP
        px("Test symmetry frakR and CP")
        for a=1:2, b=1:2, α=1:2, β=1:2, γ=1:2
            x = F[(a,b), (α,β,γ)]
            Fu = F[(ul(a),ul(b)), (α,β,γ)]

            y = conj(Fu)
            small(x-y)

            dt(q) = (-1)^ul(q)
            # ε = dt(α)*dt(β)*dt(γ)
            ε = (-1)^(α+β+γ+1)
            small(x-ε*Fu)
        end

        px("Test symmetry R_{2π/3}")
        W = get_W(;dims=3)
        for a=1:2, b=1:2
            FF = F_vecs(a,b,de;R=R)
            y = cis((a-b)*2π/3)
            A_ou_B = W-8*y*I
            M = A_ou_B*FF
            px("Norm ",norm(FF)," ",norm(A_ou_B)," ",norm(M))
        end

        px("Test symmetry R_{2π/3} on equivalent system")
        for a=1:2, b=1:2
            FF = F_vecs(a,b,de;R=R)
            y = cis((a-b)*2π/3)
            C, Ctilde = get_tilde(y)
            # if (a,b) == (1,2)
                # B, Btilde = get_tilde(conj(y))
                # Btilde = conj.(Btilde)
            # end
            # px("TRES BIZARRE QUE CE SOIT PAREIL POUR 12 ET 21 !!!")
            P = 3*2^5*Ctilde
            x = norm(C*(W-8*y*I)-P)
            small(x)
            M = P*FF
            to_tex(Ctilde;name=string("a,b=",a,",",b," C_tilde"))
            px("a,b=",a,",",b," Norm ",norm(FF)," ",norm(P)," ",norm(M)," ",x)
        end
        x = abs(F_vecs(2,1,de;R=R)[4] + im*F_vecs(2,1,de;R=R)[5])
        println("X=",x)
        x = abs(γ1 + im*F_vecs(2,1,de;R=R)[5])
        println("X=",x)
    end
end
