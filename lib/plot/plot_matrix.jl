import CairoMakie

create_dir_loc(path) = if !isdir(path) mkdir(path) end

function create_mat(M)
    B = copy(M)
    N_rows = size(M,1)
    for i=1:N_rows
        B[i,:] = M[N_rows-i+1,:]
    end
    A = transpose(B)
    abs.(A)
end

function save_fig(fig,name,dir)
    suff = dir[end]=="/" ? "" : "/"
    path = "./"*dir*suff
    create_dir_loc(path)
    full_path = path*"mat_"*name*".pdf"
    CairoMakie.save(full_path,fig)
    println("Wavefun plotted at ",full_path)
end

function plot_mat(M;name="",dir="plots")
    A = create_mat(M)
    fig = CairoMakie.heatmap(A)
    save_fig(fig,name,dir)
end

function plot_mats(M;name_dir="",dir="plots",titles="")
    n = length(M)
    fig = Figure(size = (300*n, 300))
    for i=1:n
        A = create_mat(M[i])
        ax, hm = CairoMakie.heatmap(fig[1,i], A)
        if titles != ""
            ax.title = titles[i]
        end
    end
    save_fig(fig,name_dir,dir)
end


function test()
    M = [1 2 3 4 5;
         6 7 8 9 10]
    plot_mat(M)
end

# test()
