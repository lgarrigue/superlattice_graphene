include("solvers/lobpcg.jl")
import LaTeXStrings
import CairoMakie
px = println

mutable struct BandDiagram
    # Bands diagram
    Klist # in reduced coordinates
    Knames
    n_bands # number of eigenvectors we compute
    k2lin # preconditionner for LOBPCG, square of ks, in linear list
    tol # in LOBPCG
    solver
    as_star
    n_if_action

    # Data
    Knames_to_k
    k_to_E
    ik_j_to_E
    k_to_φ
    k_to_ik
    ik_to_k
    dim

    # Parameters for plots
    root_path
    full_path
    resolution_bands
    energy_unit_plots # ∈ ["eV","Hartree"]
    path_bandwidths
    displayed_info # list of strings
    indicate_number_of_bands # if displays n_bands_window on the plot
    graph_names
    plot_legend
    plot_yticks
    numbering_bands # if displays the number of each band around it

    # Limits in y
    plot_all_bands
    energy_center
    energy_amplitude
    n_bands_window # number of bands in the window, for each graph

    function BandDiagram(;dim=2)
        bd = new()
        bd.root_path = "./plots_bands/"
        bd.full_path = ""
        # bd.folder_plots_bands = "bands"
        bd.resolution_bands = 5
        bd.tol = 1e-6
        bd.solver = "LOBPCG"
        bd.dim = dim

        K1 = [-2,1]/3; K2 = [-1,2]/3
        Γ = [0.0,0.0]; Γ2 = 2*K1-K2; M = Γ2/2
        bd.Klist = [K2,K1,Γ2,M,Γ]
        bd.Knames = ["K_2","K_1","Γ'","M","Γ"]

        bd.n_bands = 5
        bd.displayed_info = []
        bd.indicate_number_of_bands = false
        bd.graph_names = []
        bd.plot_legend = false
        bd.plot_yticks = true
        bd.n_bands_window = []
        bd.numbering_bands = true

        # Limits in y
        bd.plot_all_bands = true
        bd.energy_center = 0
        bd.energy_amplitude = 0
        bd
    end
end

function add_lattice_vectors(a1,a2,bd)
    bd.as_star = a_star_from_a_local(a1,a2)
end

create_dir_loc(path) = if !isdir(path) mkdir(path) end

######################### Core solve eigenmodes

# Applies LOBPCG
function apply_lobpcg(H,l,bd,X0,precond;maxiter=100,tol=bd.tol)
    # px("Applies LOBPCG")
    (λs,cv,Xf) = solve_lobpcg(H,l,precond;maxiter=maxiter,tol=tol,X0=X0,full_diag=true)
    (λs,cv,Xf)
end

# Solve H_K for one K
# X0 is the starting vector for LOBPCG
function solve_one(H,n,bd,precond,X0=-1,always_random=false)
    Xf = -1; E = -1
    if bd.solver=="LOBPCG"
        X = []
        if X0 == -1 || always_random
            X = randn(ComplexF64,n,bd.n_bands)
        else
            X = X0
        end
        # X = randn(ComplexF64,n,bd.n_bands)
        (E,c,Xf) = apply_lobpcg(H,bd.n_bands,bd,X,precond)
    else
        (E,φs) = eigen(H)
        E = E[1:bd.n_bands]
        # φs = φs[:,1:bd.n_bands] # Peut etre pas bon pour φs mais on s'en sert pas
        # px(E)
    end
    (E,Xf)
end

######################### Computes the band diagram

function get_eigenvalues(k_red,Hop,Xi,bd;method="k_to_matrix")
    if method=="k_to_matrix"
        Hk = Hop(k_red)
        # px("TYPE ",typeof(Hk))
        # px("K ",K0," ",K1)
        # test_hermitianity(Hk,"Hk")
        # test_part_hole_sym_matrix(Hk,bd,"Hk")
        h = norm(Hk-Hk')/norm(Hk)
        if h>1e-4 px("Be careful, H not exactly Hermitian : ",h) end
        # px("TYPE hk ",typeof(Hk))
        H = (Hk+Hk')/2
        # px("TYPE hk2 ",typeof(Hk))
        # H = H # H = Hermitian(H) CHANGES THE TYPE !!!
        # px("TYPE h ",typeof(H))
        precond = bd.k2lin(k_red)
        n = size(H,1)
        (E,X) = solve_one(H,n,bd,precond,Xi)
        # Selects eigenvalues around the Fermi level
        E = E[1:bd.n_bands]
        return (E,X)
    elseif method=="k_to_action" # if we give the action
        L = Hop(k_red)
        precond = bd.k2lin(k_red)
        (E,X) = solve_one(L,bd.n_if_action,bd,precond,Xi)
        # Selects eigenvalues around the Fermi level
        E = E[1:bd.n_bands]
        return (E,X)
    elseif method=="k_to_E"
        E = Hop(k_red)
        # px("type ",typeof(E)," ",length(E))
        n = min(length(E),bd.n_bands)
        E = E[1:n]
        return (E,-1)
    end
end

function init_klist_knames(Klist,Knames,bd)
    bd.Klist = Klist
    bd.Knames = Knames
    bd.Knames_to_k = Dict()
    for i=1:length(Klist)
        bd.Knames_to_k[Knames[i]] = Klist
    end
end

macro maybe_threads(flag, expr) # search "conditional threading julia" on google
    quote
        if $(flag)
            Threads.@threads $expr
        else
            $expr
        end
    end |> esc
end

# Computes spectrum for eigenvalues between l1 and n_bands
# It is paralellized
# H is H(q), function depending on q and giving a matrix

function init_ik_j_to_E(np,n,bd)
    if bd.ik_j_to_E == -1
        bd.ik_j_to_E = zeros(np, n)
    end
end

function spectrum_on_a_path(Hop,bd;print_progress=false,method="k_to_matrix",multithread=false)
    res = bd.resolution_bands
    n = length(bd.Klist)
    n_path_points = res*n
    X = -1
    # σ = zeros(n_path_points, bd.n_bands)
    if print_progress px("Computation of the band diagram. Progress in % (multi-threaded) :") end
    bd.ik_j_to_E = -1
    bd.k_to_ik = Dict()
    bd.ik_to_k = Dict()
    bd.k_to_E = Dict()
    bd.k_to_φ = Dict()
    for ik=1:n
        K0 = bd.Klist[ik]
        K1 = bd.Klist[mod1(ik+1,n)]
        path = [(1-t/res)*K0 .+ (t/res)*K1 for t=0:res-1]
        @maybe_threads multithread for s=1:res
        # Threads.@threads for s=1:res
        # for s=1:res
            k_red = path[s]
            # px("Here ",s)
            # for s=1:res
            (Es,X) = get_eigenvalues(k_red,Hop,X,bd;method=method)
            bd.k_to_E[k_red] = Es
            if method!="k_to_E"
                bd.k_to_φ[k_red] = X_to_φs(X)
            end
            indice = (ik-1)*res + s
            # σ[indice,:] = E
            init_ik_j_to_E(n_path_points,length(Es),bd)
            bd.ik_j_to_E[indice,:] = Es
            bd.ik_to_k[indice] = k_red
            bd.k_to_ik[k_red] = indice
            # does not work because it's parallelized
            if print_progress
                percentage = round(100*(indice-1)/(n*res), digits=0)
                print(percentage," ")
            end
        end
    end
    px("\n")
    bd.ik_j_to_E
end



# get the min and max of the band diagrams
function extrema_σs(σs)
    m = minimum([minimum(σ) for σ in σs])
    M = maximum([maximum(σ) for σ in σs])
    m, M
end

######################### Plots band diagram

function a_star_from_a_local(a1,a2) # vector in cartesian direct space to cartesian Fourier space
    lat = 2π*inv([a1 a2]')
    a1_star = lat[1:2,1]
    a2_star = lat[1:2,2]
    a1_star,a2_star
end

# vector in reduced Fourier space to cartesian Fourier space
function k_red2cart_local(k,bd)
    if bd.dim==1
        return 1
    end
    k[1]*bd.as_star[1] + k[2]*bd.as_star[2] 
end

# From the numbers of the band diagram, produces a plot of it
# σs contains a list of σ, which will be displayed in different colors. σ is a matrix, first index is the index of the path, second one is the index of the energy
function plot_band_diagram(σs,name,bd;post_name="",colors=fill(:black,100),extension="pdf",linewidth=0.5,linestyles=fill(nothing,100),resolution=(1000,700),export_in_heavy_stuff_dir=true) # zero_central_energies : if true, the central energies will be at 0 artificially
    σ0 = filter(x-> x!= nothing,σs)
    # Prepares the lists to plot
    n = length(bd.Klist)
    res = bd.resolution_bands
    n_path_points = res*n
    
    # Limits in y
    ylims = bd.energy_center - bd.energy_amplitude , bd.energy_center + bd.energy_amplitude
    if bd.plot_all_bands
        ylims = extrema_σs(σ0)
    end

    # Computes k-path
    lengths_paths = [norm(k_red2cart_local(bd.Klist[mod1(i+1,n)],bd) .- k_red2cart_local(bd.Klist[i],bd)) for i=1:n]
    lengths_paths /= sum(lengths_paths)
    if bd.dim==1
        lengths_paths = 1
    end
    dx_list = lengths_paths/res

    x_list = [] # list of x, describes each momentum
    starts_x = []; end_x = 0
    for j=1:n
        dx = dx_list[j]
        cur_l = end_x .+ [dx*(i-1) for i=1:res]
        x_list = vcat(x_list,cur_l)
        push!(starts_x,end_x)
        end_x += res*dx
    end

    # Builds figure
    f = CairoMakie.Figure(;size=resolution)
    ax = CairoMakie.Axis(f[1, 1])#, ylabel="meV?")

    start_x = x_list[1]
    end_x = x_list[end] + dx_list[end] # fictitious x point, to loop the path
    px("xlims ",start_x," ",end_x)
    CairoMakie.limits!(ax, start_x, end_x, ylims[1], ylims[2]) # x1, x2, y1, y2
    # ax.yticks = (ylims[1] : 50 : ylims[2])
    dy = (ylims[2]/ylims[1])/10
    # ax.yticks = (ylims[1]:dy:ylims[2])
    # pl = Plots.plot(size=(1000,1100),ylims=ylims,legend=false) #:topright)
    
    # Band diagram
    σss = σ0
    n_xs = length(x_list)
    n_graphs = length(σss)
    lines = []
    for g=1:n_graphs
        n_l = length(σss[g][1,:])
        for l=1:n_l
            points = [CairoMakie.Point2f(x_list[i],σss[g][i,l]) for i=1:n_path_points]
            points = vcat(points,[CairoMakie.Point2f(end_x,σss[g][1,l])])
            line = CairoMakie.lines!(points,color=colors[g],linewidth=linewidth,linestyle=linestyles[g])
            if l==1
                push!(lines,line)
            end
            # Numbering the bands
            if bd.numbering_bands 
                # Get the x the the center of the label
                jj = floor(Int,n_xs/2)
                center = x_list[jj] 
                y = σss[g][jj,l]
                # Shift
                direction = mod(g,2)==0 ? 1 : -1
                shift = direction*mod1(l,6)*end_x/100
                ds = direction*end_x/600
                # Write text
                fs = 5
                CairoMakie.text!(center+shift+ds,y; text=string(l),fontsize=fs,color=colors[g])
                # Write small lines
                points = [CairoMakie.Point2f(center,y), CairoMakie.Point2f(center+shift,y)]
                line_text = CairoMakie.lines!(points,color=colors[g],linewidth=0.5)
            end
        end
    end

    # Horizontal axis : vertical lines and ticks
    list_x_vert_labels = vcat([x_list[(i-1)*res+1] for i=1:n],[x_list[end]+dx_list[end]])
    m = length(list_x_vert_labels)
    labels = [LaTeXStrings.latexstring(string("\$",bd.Knames[mod1(i,n)],"\$")) for i=1:n+1]

    glob_shift = -0.06
    function shft(i,m)
        if i==1
            return 0.01
        elseif i==m
            return -0.08
        else
            return glob_shift
        end
    end
    pos = [CairoMakie.Point2f(list_x_vert_labels[i] + shft(i,m), ylims[1]-abs(ylims[1])*0.0) for i=1:m]
    CairoMakie.text!(labels, position=pos, fontsize=80,font ="Arial bold")
    ax.xticks = (list_x_vert_labels,["" for i=1:m])

    # Vertical axis : ticks
    ax.yticklabelsize = 40
    if !bd.plot_yticks # hide ticks
        ax.yticklabelsize = 0
    end
    ax.yticklabelpad = -90 # pull tick labels inward
    ax.yticks = round(ylims[1],digits=1):0.2:round(ylims[2],digits=1)

    # Counts mean number of bands displayed
    for g=1:n_graphs
        c = []
        M = length(σss[g][:,1])
        for m=1:M
            E = filter(x -> ylims[1] < x < ylims[2], σss[g][m,:])
            push!(c,length(E))
        end
        m = sum(c)/length(c)
        count = round(m,digits=1)
        push!(bd.n_bands_window,count)
        # px("Graph ",g," number of bands in window ",length(E))
    end

    # More information on the bottom left
    if bd.displayed_info != []
        infos = bd.displayed_info
        if bd.indicate_number_of_bands 
            s = "Nbands of "
            for g=1:n_graphs
                s=s*string(bd.graph_names[g],"=",bd.n_bands_window[g]," ")
            end
            push!(infos,s)
        end
        nl = length(infos)
        ys = [ylims[1]+(7-i)*(ylims[2]-ylims[1])/30 for i=1:nl]
        x0 = list_x_vert_labels[1]
        x1 = list_x_vert_labels[end]
        x = x0 + (x1-x0)/10
        pos = [CairoMakie.Point2f(x,y) for y in ys]
        CairoMakie.text!(infos, position=pos, fontsize=30, font ="Arial bold")
    end

    # Legend
    if bd.plot_legend
        CairoMakie.axislegend(ax, lines, bd.graph_names)
    end

    # Saves
    s = string(name,"000000000")
    title = s[1:min(6,length(s))]
    title = ""
    paths_plots = []
    # if typeof(bd.root_path)==String # case bd.root_path is a path, and there is only one
        # paths_plots = [bd.root_path]
    # else # case bd.root_path is a list of several paths
        # paths_plots = bd.root_path
    # end
    ROOT = "/home/louis/heavy_stuff/graphs_code/"
    path = "./"
    if export_in_heavy_stuff_dir
        path = ROOT
    end
    path = ROOT*bd.root_path
    if bd.full_path!=""
        path = bd.full_path
    end
    create_dir_loc(path)

    # for path in paths_plots
    # create_dir_loc(path)
    path_plot = path*"bd"*title*"_"*post_name*"."*extension
    CairoMakie.save(path_plot,f)
    px("################ Plot saved in ",path_plot)
    # end
end

######################### Plots path

function plot_sequence_of_points(Klist,Knames,bd;shifts_text=0,color=:black,linewidth=1,dots=true)
    n = length(Klist)
    Kcart = [k_red2cart_local(Klist[i],bd) for i=1:n]
    Kxs = [Kcart[i][1] for i=1:n]
    Kys = [Kcart[i][2] for i=1:n]

    # Preparation 
    points = [CairoMakie.Point2f(Kxs[i],Kys[i]) for i=1:n]

    # Lines
    loop = vcat(points,[points[1]])
    CairoMakie.lines!(loop,color=color,linewidth=linewidth)

    # Annotations
    st = shifts_text
    if shifts_text==0 st = [[0.0,0.0] for i=1:n] end
    points_annotations = [CairoMakie.Point2f(Kxs[i]+st[i][1],Kys[i]+st[i][2]) for i=1:n]
    filtered_names = []; filt_pos = []
    for i=1:n
        s = Knames[i]
        if s != ""
            push!(filtered_names,s)
            push!(filt_pos,points_annotations[i])
        end
    end
    Knames_tex = [LaTeXStrings.latexstring(string("\$",filtered_names[i],"\$")) for i=1:length(filtered_names)]
    # CairoMakie.text!(Knames_tex, position=filt_pos,textsize=30)
    for i=1:length(filtered_names)
        CairoMakie.text!(Knames_tex[i], position=filt_pos[i],fontsize=50)
    end

    # Points
    if dots
        CairoMakie.scatter!(points,color=:black)
    end
end

plot_one_vector(v,name,bd;shift_text=[0.0,0.0],linewidth=1,color=:black) = plot_sequence_of_points([[0.0,0],v],["",name],bd;shifts_text=[[0.0,0.0],shift_text],linewidth=linewidth,color=color)

# Plots the path in momentum space used for computing the band diagrams
function plot_path(Klist,Knames,bd)
    # Init
    res = 1000
    f = CairoMakie.Figure(;size=(res,res))

    ax = CairoMakie.Axis(f[1, 1],aspect=1)
    ax.aspect = CairoMakie.AxisAspect(1)

    CairoMakie.hidedecorations!(ax)
    CairoMakie.hidexdecorations!(ax, grid  = false)
    CairoMakie.hideydecorations!(ax, ticks = false)

    lims = 4
    CairoMakie.limits!(ax, -lims, lims, -lims, lims) # x1, x2, y1, y2

    # q1, q2, q3
    q1 = [-1,-1]/3; q2 = [2,-1]/3; q3 = [-1,2]/3
    plot_one_vector(q1,"q_1",bd;shift_text=[0,-0.25],linewidth=3)
    plot_one_vector(q2,"q_2",bd;shift_text=[-0.1,-0.25],linewidth=3)
    plot_one_vector(q3,"q_3",bd;shift_text=[-0.1,0.1],linewidth=3)

    # Hexagons
    hexs = []
    hexag = [q3,q3+q2,q2,q2+q1,q1,q3+q1]
    push!(hexs,hexag)
    push!(hexs,[hexag[i]+[1,0] for i=1:length(hexag)])
    push!(hexs,[hexag[i]+[0,1] for i=1:length(hexag)])
    push!(hexs,[hexag[i]+[1,0] for i=1:length(hexag)])
    push!(hexs,[hexag[i]+[-1,1] for i=1:length(hexag)])
    push!(hexs,[hexag[i]+[-1,0] for i=1:length(hexag)])
    for hex in hexs
        plot_sequence_of_points(hex,["" for i=1:length(hex)],bd;linewidth=1,dots=false,color=:grey)
    end

    # Plot list
    shifts_text = [[0.0,0.0] for i=1:length(Knames)]
    for i=1:length(Knames)
        if Knames[i]=="K_2" shifts_text[i] = [-0.5;-0.3] end
        if Knames[i]=="K_1" shifts_text[i] = [-0.5;-0.1] end
        if Knames[i]=="Γ"   shifts_text[i] = [0.1;-0.2] end
        if Knames[i]=="M"   shifts_text[i] = [-0.3;0] end
        if Knames[i]=="Γ'"  shifts_text[i] = [-0.5;-0.4] end
    end
    plot_sequence_of_points(Klist,Knames,bd;color=:blue,linewidth=5,shifts_text=shifts_text)

    # a1*,a2*
    plot_one_vector([-1.0,0.0],"a_{1,M}^*",bd;shift_text=[0,-0.3],color=:orange,linewidth=2)
    plot_one_vector([0.0,-1.0],"a_{2,M}^*",bd;shift_text=[0.1,-0.1],color=:orange,linewidth=2)

    # Saves
    path_local = "band_diagrams_bm_like/"
    paths_plots = [path_local]

    for path in paths_plots
        save(path*"path_bands_diagram.pdf",f)
    end
    px("Momentum path in dual space plotted")
end

function spectrum_on_a_path_1d(k2E_map,bd;print_progress=false,method="k_to_matrix",multithread=false)
    res = bd.resolution_bands
    X = -1
    kreds = [j/res for j=0:res] # red and not cart !
    Nk = length(kreds)
    for ik=1:Nk
        kred = kreds[ik]
        Es = k2E_map(kred)
        n_l = Int(length(Es))
        if !isnan(bd.n_bands) n_l = bd.n_bands end
        if ik==1 bd.ik_j_to_E = zeros(Nk, n_l) end
        bd.ik_j_to_E[ik,:] = Es[1:n_l]
    end
    bd.ik_j_to_E
end

function plot_band_diagram_1d(σss,name,bd;post_name="",colors=fill(:black,100),numbering_bands=false,extension="pdf",linewidth=0.5,linestyles=fill(nothing,100),resolution=(1000,700),center=0,amplitude=100,ground=false,export_in_heavy_stuff_dir=true) # zero_central_energies : if true, the central energies will be at 0 artificially
    colors = [:red,:blue,:black]
    res = bd.resolution_bands
    kreds = [j/res for j=0:res]
    Nk = length(kreds)
    # Limits in y
    # ylims = bd.energy_center - bd.energy_amplitude , bd.energy_center + bd.energy_amplitude
    f = CairoMakie.Figure(;size=resolution)
    ax = CairoMakie.Axis(f[1, 1])#, ylabel="meV?")
    # CairoMakie.limits!(ax, 0, 1, ylims[1], ylims[2]) # x1, x2, y1, y2
    ymax = center+amplitude
    ymin = center-amplitude
    # if ground
        # ymin = minimum(minimum(Ess))
    # end
    CairoMakie.limits!(ax, 0, 1.1, ymin, ymax) # x1, x2, y1, y2
    
    # Band diagram
    n_graphs = length(σss)
    lines = []
    for g=1:n_graphs
        n_l = length(σss[g][1,:])
        # n_bands = minimum(bd.n_bands,n_l)
        for l=1:n_l
            points = [CairoMakie.Point2f(kreds[i],σss[g][i,l]) for i=1:Nk]
            line = CairoMakie.lines!(points,color=colors[g])#,linewidth=linewidth,linestyle=linestyles[g])
        end
    end

    # Saves
    path = "./"
    if export_in_heavy_stuff_dir
        path = "../../../heavy_stuff/graphs_code/improved_graphene/"
        create_dir_loc(path)
    end
    path_plot = path*"bd_"*post_name*"."*extension
    CairoMakie.save(path_plot,f)
    px("Plot saved in ",path_plot)
end

function plot_spectrum(Ess_C;center=0,amplitude=1,ground=false,draw_lines=false)
    N = length(Ess_C)
    if norm(imag.(Ess_C))>1e-4
        px("SPECTRUM NOT REAL")
    end
    Ess = [real.(Es) for Es in Ess_C]
    res = 1000
    f = CairoMakie.Figure(;size=(500,res))
    ax = CairoMakie.Axis(f[1, 1])
    ymax = center+amplitude
    ymin = center-amplitude
    if ground
        ymin = minimum(minimum(Ess))
    end
    dx = 0.2 # separation in x between two spectrums
    CairoMakie.limits!(ax, -2, dx*(N+1), ymin, ymax) # x1, x2, y1, y2
    # Points if not lines
    points = []
    for p=1:N
        points = vcat(points,[CairoMakie.Point2f(p*dx,Ess[p][i]) for i=1:length(Ess[p])])
    end
    CairoMakie.scatter!(points,color=:black,markersize=3)
    # Lines if lines
    M = length(Ess[1]) # all are considered to have the same length
    lines = []
    for m=1:M
        points = [CairoMakie.Point2f((i-1)*dx,Ess[i][m]) for i=1:N]
        push!(points,CairoMakie.Point2f(N*dx,Ess[1][m]))
        if draw_lines
            line = CairoMakie.lines!(points)
        end
    end
    # Add numerotation of spectrum
    for j=1:length(Ess[1])
        CairoMakie.text!(-0.5,Ess[1][j];text=string(j),fontsize=10,color=:black)
    end
    # Saves
    path = "./"
    # create_dir_loc(path)
    path_plot = "spectrum.pdf"
    CairoMakie.save(path_plot,f)
    px("Plot of spectrum saved in ",path_plot)
end
