include("../0_low_level/operations_on_fc_functions.jl")
include("../1_common_graphene_manips/common_graphene.jl")

include("../lib/band_diagram.jl")
include("../lib/pt/asymptotic_behavior.jl")
include("../lib/gram_schmidt.jl")

include("../2_monolayer_derivatives/graphene_derivatives.jl")

include("../3_effective_system/compute_i_dmode.jl")
include("../3_effective_system/stock_module_coefs.jl")
include("../3_effective_system/stock_module_solutions.jl")
include("../3_effective_system/solve_exact.jl")
include("../3_effective_system/solve_effective.jl")
include("../3_effective_system/plot_curves.jl")
include("../3_effective_system/misc.jl")

include("../lib/stocks/stock_module.jl")

mean(l) = sum(l)/length(l)


M_DIR, K_DIR = get_MK()

function continuously_get_next_eigenvector(last_ψ,var,oc,St_sol;previous_is=[])
    i = -1
    if length(previous_is) >= 1
        i = previous_is[end]
    else
        px("COMPUTES FULL I DMODE")
        i = get_i_dmode(oc,St_sol)
    end
    px("i_dmode is ",i)
    ψs, is = get_ψs_around(i,oc,St_sol; amplitude = 200)
    ψ = get_ψ(i,oc,St_sol) # to define it globally
    # Searches for the right ψ
    if haskey(last_ψ,oc.e) && !(var in ["invep"])
        lψ = last_ψ[oc.e]
        dists = [dist_quo(φ,lψ) for φ in ψs]
        ds = []
        for a=1:length(is)
            if dists[a] < 1
                push!(ds,(is[a],dists[a]))
            end
        end
        px("distances are ",ds)
        ii = argmin(dists)
        px("min dist is ",dists[ii]," for ",ii)
        ψ, i = ψs[ii], is[ii]
    end
    px("new_i_dmode is ",i)
    ψ, i
end

function produce_line(λs,var,oc,St_sol; V_base=0)
    l = []; is_ex = []; is_ef = []; new_λs = []
    last_ψ = Dict()
    last_dist = -1
    for j=1:length(λs)
        λ = λs[j]
        oc = modify_oc(λ,var,oc; V_base=V_base)
        px("Computation at invε=",oc.invε,", kred=",oc.kred,", V=",oc.V_intensity)
        ##### Exact
        # first do exact because for effective, some numbers are reused from the exact computation
        oc.e = "ex"
        px("EXACT IDMODE")
        ψ_ex, i_ex = continuously_get_next_eigenvector(last_ψ,var,oc,St_sol;previous_is=is_ex)
        ##### Effective
        oc.e = "ef"
        px("EFFECTIVE IDMODE")
        ψ_ef, i_ef = continuously_get_next_eigenvector(last_ψ,var,oc,St_sol;previous_is=is_ef)
        px("i EX ",i_ex,"i EF ",i_ef)
        # px("i_ef is ",i_ef," ",oc.order,"   -   ",oc.k_dep)
        # px("i_ex is ",i_ex," ",oc.order,"   -   ",oc.k_dep)

        ##### Compare
        # dists = [dist_quo(ψ_ex,ψ_ef) for ψ_ef in ψs_ef] # compares with also "excited" states because it can happen that a higher level comes to a lower and makes a band crossing
        # d = dist_quo(ψ_ex,ψ_ef)
        # i_ef = argmin(dists)
        # i = get_i_dmode(oc,St_sol)-3
        # α = get_α_ef(i,oc,St_sol)
        # norms = [norm(α[j]) for j=1:length(α)]
        # px("\n##### ",norms,"\n")
        # px("\n##### ",norms*oc.invε,"\n")
        # for j in [1]
            # d = abs(get_α_ef(i,oc,St_sol)[j])
            # px("\n---------> For ε=",1/oc.invε," |α|",j,"=",d,"\n")
        # end
        # for j in [3]
            # d = abs(get_α_ef(i,oc,St_sol)[j])*oc.invε # <---------- REMOVE THIS, JUST FOR AN EXPERIMENT
            # px("\n---------> For ε=",1/oc.invε," |α|",j,"=",d,"\n")
        # end
        px("IL FAUT PLOTTER TOUTE LA NORM DE α")
        d = dist_quo(ψ_ef,ψ_ex)
        if j==1
            last_dist = d
        end
        RATIO_LIMIT = λ < 7 ? 10 : 3
        if d/last_dist < RATIO_LIMIT # if not, there is a singularity, i.e. we did not get the right eigenvector, we have a mixed one due to eigenvalue crossing
            push!(l,d); push!(is_ex,i_ex); push!(is_ef,i_ef); push!(new_λs,λ)
            last_ψ["ef"] = ψ_ef; last_ψ["ex"] = ψ_ex
            last_dist = d
        end

        # px("---> For λ=",λ," i_ex=",i_ex," i_ef=",i_ef)
    end
    # Release
    new_λs, l, is_ex, is_ef
end

abscisse = Dict()
abscisse["N_macro"] = "Nmacro"
abscisse["lambda_k"] = "\\mu"
abscisse["V"] = "\\lambda"
abscisse["invep"] = "\\varepsilon"
abscisse["V_and_ep"] = "V_ep"

function one_study()
    var = "N_macro" # ∈ ["N_macro","lambda_k","V"]
    var = "lambda_k"
    var = "V"
    # var = "invep"
    λs = []
    λkred = 1e-7
    # λkred = 0.3
    short = false
    finalized = true
    oc = OneComp()
    oc.precision = short ? "normal" : "precise" # if "normal", precisions cannot go below 1e-2
    oc.invε = 7 # DONT CHOOSE 4 OR OTHER VALUES WHERE THERE IS DEGENERACY ! For 5 there is a problem with k-dep and ℓ=2 so let is like that
    if var=="V"
        if short
            res = 3
            oc.invε = 2
        end
        res = 15
        λs = get_λs(-4,-1.1;res=6) # V varies
        λs = get_λs_log_lin(-1,0.9;res=res) # V varies
        sp = []
        λs = vcat(sp,λs)
        oc.N_macro = 15
        # oc.precision = "ultra"
    elseif var=="lambda_k"
        res = 4
        if finalized
            res = 5
        end
        λs = get_λs_log_lin(0.5,1;res=5) # lambda_k
        # λs = get_λs(0.5,1.1;res=res) # lambda_k
        sp = [1e-2,0.05,0.1,0.2,0.5,0.7,0.8,1,1.5,1.8,2,2.5,8.7,9.4,9.7,9.9,1.02]
        λs2 = get_λs_log_lin(1.03,1.15;res=5) # lambda_k
        λs = vcat(sp,λs)
        λs = vcat(λs2,λs)
        oc.N_macro = 15
        # λs = sp
    elseif var in ["invep"]
        λs = [8,11,14,17,20,23] # ε
        λs = [7,9,11] # ε
        λs = [19,16,13,10,7,4] # ε
        λs = [20,17,14,11,8,5] # ε
        λs = [1,2,3,4,5,6,7,8]
        λkred = 0.5
        if short
            λs = [1,2,3,4,5,6,7]
        end
    end
    sort!(λs)
    px("λs \n",λs)

    px("STUDY "*var)
    St_sol = create_stock_sol()
    # λkred = 0.5
    oc.kred = λkred*K_DIR
    oc.V_intensity = 0
    oc.n_zeroth = 2
    # oc.N_macro = 15
    # oc.order = 1
    # oc.k_dep = "k_independent"
    V_base = 0 # for var="V_and_ep" only

    # Build oc list
    ocs = []
    orders = [0,1,2]
    if !finalized
        orders = [0]
    end
    # orders = [0]
    for order in orders, kdep in (var=="V" ? kdo(order) : ["k_independent"])
    # for order in orders, kdep in kdo(order)
        oc0 = change_oc(order,"order",oc)
        oc1 = change_oc(kdep,"k_dep",oc0)
        push!(ocs,oc1)
    end
    for zs in [6]
        oc3 = change_oc(0,"order",oc)
        oc4 = change_oc("k_independent","k_dep",oc3)
        oc5 = change_oc(zs,"n_zeroth",oc4)
        if finalized
            push!(ocs,oc5)
        end
    end

    filename = var*"_varies"
    # Legend
    variable_string = LaTeXString(abscisse[var])
    leg_info = []
    # Plots
    miny = finalized ? 1e-4 : 1e-7

    plot_curves(λs,var,ocs,St_sol,filename,variable_string,leg_info;
                miny=miny,V_base=V_base,scat=var=="invep",orientation=:vertical,finalized=finalized)
end

one_study()
