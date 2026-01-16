using CairoMakie
include("../0_low_level/my_solver_graphene.jl")
include("../1_common_graphene_manips/common_graphene.jl")

function test_equality()
    # DFTK stuff
    # N = 32; Nz = 192; gauge_param = 1
    p = get_standard_param_module("dftk_import")
    u1_fc = u_fcl_n(1,p)
    b = p.basis

    resi = apply_schro(u1_fc,p.v_dir,b;kred=b.K_red_2d) - E_fermi(p)*u1_fc
    for ik=1:b.N3d
        x = abs(resi[ik])
        k = ik2kred_3d(ik,b)
        kc = ik2kcart_3d(ik,b)
        akz = abs(k[3])
        nkc = norm(kc)
        if x > 0.001 && akz < 3
            # print(k," ",nkc,"-")
            a=1
            # print("kz=",k[3]," x=",x," ")
        end
    end
    # proportions_fourier_modes(resi,3,b;only_kz=true)
    # px("re ",norm(real.(resi))," im ",norm(imag.(resi)))
    # px("res ",norm(resi)/norm(u1_fc))

    vd = pot_fc(p)
    proportions_fourier_modes(vd,3,b;only_kz=false)

    proportions_fourier_modes(vd,3,b;only_kz=false)
    vd = myifft(vd,b.Vol)
    proportions_fourier_modes(resi,3,b;only_kz=false)

    n_eigenvectors = 20
    (λs,us) = diag_monolayer_at_k_mine(b.K_red_2d,vd,n_eigenvectors,b)

    # Exact
    px("Exact ",p.Es_K[1:n_eigenvectors] .- p.Es_K[6])
    # Effective
    px("Effective ",λs .- λs[5])
    
    xs = []
    ys = []
    x_ex = p.Es_K[1:n_eigenvectors].-p.Es_K[1]
    x_eff = λs.-λs[4]

    y_eff = 1*ones(length(x_eff))
    y_ex = 0*ones(length(x_ex))

    f = Figure()
    ax = Axis(f[1, 1])
    colors = [:red,:blue] # eff en blue
    xs = [x_ex,x_eff]
    ys = [y_ex,y_eff]
    for j=1:2
        for i=1:length(xs[j])
            x1 = xs[j][i]
            x2 = xs[j][i]
            y1 = ys[j][i]
            y2 = 0.5
            lines!(ax, [x1,x2], [y1,y2], color = colors[j], linewidth = 0.5)
        end
    end
    CairoMakie.save("plots/energies.pdf",f)
end

test_equality()
