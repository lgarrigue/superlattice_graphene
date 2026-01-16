using CairoMakie
include("../../lib/stocks/stock_module.jl")
px = println

# coords = (n,q), n is the n^th state we have in the bucket, and q is an extra coord, possibly just 1 if the states are k-idependent, and possibly the direction of k for instance if the states depend on the direction of k
function create_stock_micro_coefs(ms, p; gram=true)

    # 2) Define default parameters values
    pS = Dict()
    pS["gram"] = gram
    pS["ms"] = ms

    # 3) Define how to compute the quantities
    function compute_ξ_family(coords,St::Stock)
        Ns = length(ms.states)
        q = coords
        states = []
        for j=1:Ns
            sta = ms.states[j]
            ψ = get_state(sta,q)
            push!(states,ψ)
        end
        if St.p["gram"]
            return gram_schmidt_family(states;force=false)
        end
        states
        # px("Computed ξ(n)=ξ(",n,") = ",x)
    end

    function compute_ξ(coords,St::Stock)
        (n,q) = coords
        X("ξ_fam",q,St)[n]
        # sta = ms.states[n]
        # get_state(sta,q)
        # px("Computed ξ(n)=ξ(",n,") = ",x)
    end

    function compute_HmE_ξ(coords,St::Stock)
        (n,q) = coords
        ξ = X("ξ",coords,St)
        # φ = diff_op_K(ξ,"H_K-E",p)
        φ = apply_schro(ξ,p.v_dir,p.basis;kred=p.basis.K_red_3d) - E_fermi(p)*ξ
        # px("Computed hξ(n)=hξ(",φ,")")
        φ
    end

    function compute_mi∇ξ(coords,St::Stock) # (-i∇)ξ
        ξ = X("ξ",coords,St)
        φ = diff_op_K(ξ,"∇_K",p)
        # (φ1,φ2,φ3) = [-im*φ[j] for j=1:3]
        -im*φ[1], -im*φ[2]
    end

    function compute_Δξ(coords,St::Stock) # Δξ
        ξ = X("ξ",coords,St)
        diff_op_K(ξ,"Δ_K",p)
    end

    function compute_ξξ(coords,St::Stock)
        (n,m,q) = coords
        ξ = X("ξ",(n,q),St)
        φ = X("ξ",(m,q),St)
        dot(ξ,φ)
    end

    function compute_ξhξ(coords,St::Stock)
        (n,m,q) = coords
        ξ = X("ξ",(n,q),St)
        φ = X("(h-E)ξ",(m,q),St)
        s = dot(ξ,φ)
        # px(s)
        s
    end

    function compute_ξmi∇ξ(coords,St::Stock) # ξ⋅(-i∇)ξ
        (n,m,q) = coords
        ξ = X("ξ",(n,q),St)
        φ1, φ2 = X("(-i∇)ξ",(m,q),St)
        dot(ξ,φ1), dot(ξ,φ2)
    end

    function compute_Minv(coords,St::Stock)
        q = coords
        N = length(ms.states)
        M = [X("ξ(h-E)ξ", (m1+2,m2+2,q), St) for m1=1:N-2, m2=1:N-2]
        # assert_herm(M,1e-10) # ok
        inv(M)
    end

    function compute_T(coords,St::Stock)
        (a,q) = coords
        N = length(ms.states)
        T = zeros(ComplexF64,2,N-2)
        for m1=1:2, m2=1:N-2
            V = X("ξ(-i∇)ξ",(m1,m2+2,q), St)
            T[m1,m2] = V[a]
        end
        T
    end

    function compute_TMT(coords,St::Stock)
        (a,b,q) = coords
        Ta = X("T", (a,q), St)
        Tb = X("T", (b,q), St)
        px("size ",size(Ta))
        Minv = X("Minv", q, St)
        op = Ta*Minv*(Tb')
        op
    end

    # 4) Gather computation functions, same order than kinds list
    # kinds = ["ξ_fam","ξ","(h-E)ξ","(-i∇)ξ","Δξ","ξξ","ξ(h-E)ξ","ξ(-i∇)ξ","Minv","T","TMT"]
    kcf = [("ξ_fam",compute_ξ_family),
           ("ξ",compute_ξ),
           ("(h-E)ξ",compute_HmE_ξ),
           ("(-i∇)ξ",compute_mi∇ξ),
           ("Δξ",compute_Δξ),
           ("ξξ",compute_ξξ),
           ("ξ(h-E)ξ",compute_ξhξ),
           ("ξ(-i∇)ξ",compute_ξmi∇ξ),
           ("Minv",compute_Minv),
           ("T",compute_T,),
           ("TMT",compute_TMT)]
    
    # 5) Create stock
    St = Stock(kcf,pS)
    St
end

function test_plot()
    St = create_stock()
    St.p["y"] = 0.52
    # y = 0.59 # transition

    # Computes
    n_max = 100
    ns = range(1, n_max)
    zs = [X("z",n,St) for n in ns]

    # Builds a plot
    f = Figure()
    ax = Axis(f[1, 1])
    scatter!(ax, ns, zs)
    return f
end

# test_plot()
