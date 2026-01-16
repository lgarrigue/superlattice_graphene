function plot_mean_V(p) # rapid plot of V averaged over XY
    Vz = real.([sum(p.v_dir[:,:,z]) for z=1:p.basis.Nz])/p.basis.N^2
    pl = Plots.plot(Vz)
    savefig(pl,p.path_plots*"V_averaged_over_xy.png")
end

# array of Fourier coefs to cartesian direct function
function arr2fun(u_fc,p;bloch_trsf=true,k_p=-1) # u in the Fourier ball, to function in cartesian
    K_cart = kred2kcart_2d(p.basis.K_red_2d,p)
    plus_k = bloch_trsf ? K_cart : [0,0] # EQUIVALENT TO APPLY e^{iKx} !
    f(x,y,z) = 0
    for imx=1:p.basis.N, imy=1:p.basis.N, imz=1:p.basis.Nz
        mx,my,mz = p.kred_x[imx],p.kred_x[imy],p.kred_z[imz]
        if norm([mx,my,mz]) < p.plots_cutoff
            ma = mx*p.a1_star + my*p.a2_star + (k_p==-1 ? plus_k : k_p)
            g(x,y,z) = u_fc[imx,imy,imz]*cis([ma[1],ma[2],0]⋅[x,y,z])
            f = f+g
        end
    end
    f
end

function eval_f(g,res,p) # evaluates the function g
    # a = SharedArray{ComplexF64}(res,res,p.basis.Nz)
    a = zeros(ComplexF64,res,res,p.basis.Nz)
    indices = [(i,j,kiz) for i=0:res-1, j=0:res-1, kiz=1:p.basis.Nz]
    aa = p.a
    kt = p.kred_z
    # for l=1:length(indices) # 9.5s
    Threads.@threads for l=1:length(indices) # 5s
        (i,j,kiz) = indices[l]
        a[i+1,j+1,kiz] = g(i*aa/res,j*aa/res,kt[kiz])
    end
    a
end

eval_f(g,res,p) = [g(i*p.a/res,j*p.a/res,p.kred_z[kiz]) for i=0:res-1, j=0:res-1, kiz=1:p.basis.Nz]
scale_fun3d(f,λ) = (x,y,z) -> f(λ*x,λ*y,λ*z)

function simple_plot(u,fun,Z,p;n_motifs=3,bloch_trsf=true,res=25,k_p=-1)
    f = arr2fun(u,p;bloch_trsf=bloch_trsf,k_p=k_p)
    g = scale_fun3d(f,n_motifs)
    a = fun.(eval_f(g,res,p))
    b = intZ(a,p)
    # b = a[:,:,Z]
    Plots.heatmap(b,size=(1000,1000), aspect_ratio=:equal)
end

# plots a function in both direct and Fourier space, with absolute value, real and imaginary parts, at some Z coordinate
function rapid_plot(u,p;n_motifs=5,name="rapidplot",bloch_trsf=true,k_p=-1,res=25)
    Z = 5
    funs = [abs,real,imag]
    hm = [simple_plot(u,fun,Z,p;n_motifs=n_motifs,bloch_trsf=bloch_trsf,res=res,k_p=k_p) for fun in funs]
    plot_z = Plots.plot(intXY(real.(myifft(u,p.basis.L)),p))
    push!(hm,plot_z)
    size = 600
    r = length(hm)
    pl = Plots.plot(hm...,layout=(r,1),size=(size+100,r*size),legend=false)
    # pl = Plots.plot(hm...,layout=(1,r),size=(r*size,size-200),legend=false)
    full_name = string(p.path_plots,name,"_N",p.basis.N,"_Nz",p.basis.Nz)
    savefig(pl,full_name*".png")

    # Fourier
    uXY = intZ(abs.(u),p)
    uZ = intXY(abs.(u),p)
    plXY = Plots.heatmap(uXY,size=(1500,1000))
    plZ = Plots.plot(uZ,size=(1500,1000))
    pl = Plots.plot(plXY,plZ)
    savefig(pl,full_name*"_fourier.png")
    px("Plot of ",name," done")
end
