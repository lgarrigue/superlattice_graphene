COLORS = ["red","brown","blue","black","brown","green","orange","pink","cyan","grey",RGBf(0.1, 0.9, 0.4),RGBf(0.4, 0.2, 0.5),RGBf(0.7, 0.1, 0.3),RGBf(0.7, 0.54, 0.9),RGBf(0.1, 0.4, 0.5),RGBf(0.9, 0.1, 0.4)]


PATH_PLOTS = "../../../tentatives/1-superlattice_graphene/pics/curves"

function plot_curves(λs,var,ocs,St_sol,name,var_str,leg_info; miny="",orientation=:horizontal, plot=true, V_base=0, scat=false, finalized=true)
    #### Builds a plot
    size_txt = 30
    # Structure
    f = Figure(resolution=(600, 600), fontsize=size_txt)
    xlabel = var_str
    # ylabel = label(Eϕ,ℓ_str,μ_str,ν_str,"rbm")
    ax = Axis(f[1, 1], 
              xscale=Makie.log10, yscale=Makie.log10, xlabel=xlabel)#, yminorticks = [10^(float(n)) for n=-15:15], yminorticksvisible = true, yminorgridvisible = true)#, xticks = [10^(-4)])ylabel=ylabel, 
    Makie.xlims!(ax, 1e-1, maximum(λs))
    Makie.xlims!(ax, minimum(λs), maximum(λs))

    ax.xticks = LogTicks(-12:1)
    ax.yticks = LogTicks(-12:1)
    # hidedecorations!(ax, grid = false) # hides x ticks

    # Legend preparation
    lines_leg = []
    labels_leg = []

    lines_ = []; lines_data = []; λss = []
    # Fills with lines
    mminy = 1
    maxy = 2
    i = 1
    is_ex = []
    fontsize_labels_l = 8
    for oc in ocs
        linestyle = nothing
        if oc.k_dep=="k_dependent"
            linestyle = :dot
        end
        px("------------------ NEW CURVE COMPUTATION, ℓ=",oc.order," kdep=",oc.k_dep," es=",oc.n_zeroth)
        λs_loc, line, is_ex, is_ef = produce_line(λs,var,oc,St_sol; V_base=V_base)
        px("*********************** OC order ",oc.order)
        push!(λss,λs_loc)
        push!(lines_data,line)
        mminy = min(mminy,minimum(line))
        maxy = max(maxy,maximum(line))
        # wrk(x) = x <1e-14 ? 1e200 : x
        # line = [wrk(line[j]) for j=1:length(line)]
        ic = oc.order+1
        if oc.n_zeroth==6
            ic = 4
        elseif oc.n_zeroth==18
            ic = 11
        end
        color = COLORS[ic]
        if !scat
            # linewidth=approx_to_linewidth[C.approx]) # line
            l = lines!(ax, λs_loc, line, color=color, linestyle=linestyle, linewidth=3)
        else
            l = CairoMakie.scatter!(ax, λs_loc, line, color=color, markersize = 20)
        end
        push!(lines_,l)
        push!(lines_leg,l)
        suff = ""
        lab = LaTeXString("\$ \\ell=$(oc.order)\$")
        if oc.order ≥ 0

            lab = oc.k_dep=="k_dependent" ? LaTeXString("\$ \\ell=$(oc.order), k_{\\text{dep}}\$") : LaTeXString("\$ \\ell=$(oc.order), k_{\\text{indep}}\$")
            lab = oc.k_dep=="k_dependent" ? LaTeXString("\$ \\mathcal{F}_{$(oc.order),k}\$") : LaTeXString("\$ \\mathcal{F}_{$(oc.order)}\$")
        end
        if oc.n_zeroth≥3
            lab = LaTeXString(string(oc.n_zeroth)*" eigenstates states at K")
            lab = LaTeXString("\$ \\mathcal{F}^{$(oc.n_zeroth)}\$")
        end
        push!(labels_leg,lab)
        i += 1
        # Indices of eigenstates for effective
        if !finalized
            for j=1:length(λs_loc)
                st = string(is_ef[j])
                tx = CairoMakie.text!(ax, λs_loc[j], line[j]; text=st,fontsize=fontsize_labels_l)
            end
        end
    end
    used_min = miny=="" ? mminy : miny, max(maxy,1)
    # Indices of eigenstates for exact
    # if !finalized
        # for j=1:length(λs)
            # CairoMakie.text!(λs[j],1; text=string(is_ex[j]),fontsize=fontsize_labels_l,color=:blue)
        # end
    # end
    # Print asymptotic behaviors
    for j=1:length(ocs)
        oc = ocs[j]
        px(var_str)
        print_results(λss[j],lines_data[j],1,1;title=string("study ℓ=",oc.order,", kdep=",oc.k_dep,", n_zeroth=",oc.n_zeroth),comparison="last")
    end
    # Slopes
    # if !finalized
    if true
        for slope in (1:4)
            colors_slopes = ["grey","orange","cyan","pink"]
            line = [(λ/10)^slope for λ in λs]
            l = lines!(ax, λs, line, color="grey", linestyle= :dashdot, linewidth=2)
            # text!(ax, 0.1, 0.5-10^(-0.5*slope), text = LaTeXString("Slope \$ $(slope) \$"), color=colors_slopes[slope])
            push!(lines_,l)
        end
    end
    if !plot
        return 0
    end
    Makie.ylims!(ax, used_min)

    # Saves
    namefile = PATH_PLOTS*"/"*name*".pdf"
    # px(namefile)
    save(namefile,f)
    px("====> Figure saved in ",namefile)

    # Legend
    res = orientation==:horizontal ? (2000, 80) : (300, 400)
    fig_leg = Figure(resolution=res, fontsize=20)
    Legend(fig_leg[1,1], lines_leg, labels_leg, halign = :left, valign = :center, orientation = orientation)
    namefile_leg = PATH_PLOTS*"/"*name*"_legend.pdf"
    save(namefile_leg,fig_leg)
    px("====> Legend saved in ",namefile_leg)
end
