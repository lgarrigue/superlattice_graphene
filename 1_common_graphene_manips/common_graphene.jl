include("../0_low_level/common_functions.jl")
include("../0_low_level/basis_graphene.jl")
include("../0_low_level/operations_on_fc_functions.jl")
include("../0_low_level/my_solver_graphene.jl")

include("../1_graphene_dftk/graphene_dftk.jl")
include("../1_graphene_dftk/operations_on_fb_functions_dftk.jl")
include("../1_graphene_dftk/import_export_funs.jl")
include("../1_graphene_mysolver/my_2d_structs.jl")

mutable struct Params
    basis
    i_state # index of first valence state
    Es_K; us_K # all eigenstates at K, in Fourier ball
    v_dir # in direct space and cube, not linear
    type # ∈ ["dftk","mine"]
    spec

    function Params(dim)
        p = new()
        dim = 3
        p.basis = Basis(dim)
        p
    end
end

function diff_op_K(ψ,kind,p)
    φ = []
    if p.type[1:4]=="dftk"
        if kind=="∇_K"
            φ = ∇_K_dftk(ψ,p)
        elseif kind=="Δ_K"
            φ = Δ_K_dftk(ψ,p)
        elseif kind=="H_K-E"
            φ = H_fb_dftk(ψ,p) - E_fermi(p)*ψ
        end
    else
        K = p.basis.K_red_3d
        if kind=="∇_K"
            φ = ∇_k_mine(ψ,p.basis;kred=K)
        elseif kind=="Δ_K"
            φ = Δ_k_mine(ψ,p.basis;kred=K)
        elseif kind=="H_K-E"
            φ = apply_schro(ψ,p.v_dir,p.basis;kred=K) - E_fermi(p)*ψ
        end
    end
    φ
end

# type ∈ ["dftk_import","dftk_new","mine"]
function get_standard_param_module(type;N=5,V_intensity=10)
    p = ""
    if type=="dftk_import"
        # p = import_u1_u2_V(32,720)
        p = import_u1_u2_V(27,576)
        # p = import_u1_u2_V(32,720)
        @assert p.i_state == 4
    elseif type=="dftk_new"
        dim = 3
        p = Params(dim)
        p.i_state = 4
        p.spec = InfoDFTK()
    elseif type=="mine_2d"
        p = generate_some_2d_param(N,V_intensity)
    elseif type=="mine_3d"
        dim = 3
        p = Params(dim)
        p.i_state = 1
    else
        px("PROBLEM IN GET STANDARD PARAM MODULE")
        @assert false
    end
    p.type = type
    p
end

# (E-H)^(-1)
function apply_pseudo_inverse(ψ,p)
    Ef = E_fermi(p)
    # px("Starts building pseudo-inverse")
    c = 0*p.us_K[1]
    if length(p.Es_K) != p.basis.Nfull
        px("BE CAREFUL, WHEN APPLYING PSEUDO INVERSE, MAYBE NOT FULL")
    end
    for i=1:length(p.Es_K)
        if !(i in [p.i_state, p.i_state+1])
            iv = (1/(Ef-p.Es_K[i]))
            dt = dot(p.us_K[i],ψ)
            φ = iv*dt*p.us_K[i]
            c += φ
            if isnan(norm(φ)) # necessary to stabilize the computation... why ? don't know
                px("pseudo inv coef NaN ",iv," ",dt)
            end
            if norm(φ) > 1e10 # necessary to stabilize the computation... why ? don't know
                px("pseudo inv ",iv," ",dt)
            end
        end
    end
    if isnan(norm(c))
        px("pseudo inverse gives globally NaN")
        @assert false
    end
    b = norm(c) > 1e10
    if b
        px("pseudo inverse gives too high norm")
        @assert false # bug, sometimes c=NaN...
    end
    c
    # px("Ends building pseudo-inverse")
    # px("Number of states taken into account in pseudo inverse: ",length(p.Es_K)," is it enough ?")
    # px("Max energy distance ",distance(p.Es_K[1], p.Es_K[length(p.Es_K)]))
end

E_fermi(p) = p.Es_K[p.i_state]
u_K(a,p) = p.us_K[p.i_state+a-1]
N_basis(p) = length(p.us_K[1])
function u_lin2fc(ψ,p)
    if p.type[1:4]=="dftk"
        return u_fb2fc_dftk(ψ,p)
    end
    myreshape(ψ,p.basis)
end

extract_eigenmodes(Es,us,p) = Es[p.i_state], Es[p.i_state+1], us[p.i_state], us[p.i_state+1]

######################### Obtain the good orthonormal basis for the periodic Bloch functions at Dirac points u1 and u2

function R_2pi3(u,p)
    if p.type[1:4]=="dftk"
        return R_fb(u,p)
    end
    R_fc_2d(u,p.basis)
end

function M_change(u,p)
    if p.type[1:4]=="dftk"
        return M_fb(u,p)
    end
    M_fc_2d(u,p.basis)
end

function P_change(u,p)
    if p.type[1:4]=="dftk"
        return P_fb(u,p)
    end
    P_fc_2d(u,p.basis)
end

function τ_k(u,kred,p)
    if p.type[1:4]=="dftk"
        # MODIFIE, AVANT ETAIT τ_fb(u,p), VERIFIER QUE C EST BON
        return τ_fb(u,kred,p)
    end
    τ_fc_2d(u,kred,p.basis)
end

function ∇_K_common(u,p)
    if p.type[1:4]=="dftk"
        return ∇_K_dftk(u,p)
    end
    ∇_k_mine(u,p.basis;kred=p.basis.K_red_2d)
end

# Turns the degenerate eigenvectors u1 and u2 in SU(2) so that the physical associated eigenvectors respect the symmetries R ϕj = ω^j ϕj. See documentation to see the details on this operation
function get_natural_basis(u1,u2,p)
    τau = cis(2π/3)
    (Ru1,Tu1) = (R_2pi3(u1,p), τ_k(u1,shift_K_R(p.basis),p))
    (Ru2,Tu2) = (R_2pi3(u2,p), τ_k(u2,shift_K_R(p.basis),p))
    d1 = Ru1.-τau*Tu1
    d2 = Ru2.-τau*Tu2

    c = (norm(d1))^2
    s = d1'*d2
    f = (c/abs(s))^2

    U = (s/abs(s))/(sqrt(1+f))*u1 + (1/sqrt(1+1/f))*u2
    V = (s/abs(s))/(sqrt(1+f))*u1 - (1/sqrt(1+1/f))*u2

    (RU,TU) = (R_2pi3(U,p),τ_k(U,shift_K_R(p.basis),p))
    (RV,TV) = (R_2pi3(V,p),τ_k(V,shift_K_R(p.basis),p))

    I = argmin([norm(RU.-τau *TU),norm(RV.-τau *TV)])
    u1_new = I == 1 ? U : V
    # px("noRm ",norm(u1))

    # ϕ2(x) = conj(ϕ1(-x)) ⟹ u2(x) = conj(u1(-x))
    # conj ∘ parity is conj in Fourier
    u2_new = conj.(u1_new)
    return u1_new, u2_new 
end

function get_natural_basis_u1_and_u2(p)
    u1, u2 = u_K(1,p), u_K(2,p)
    u1n, u2n = get_natural_basis(u1,u2,p)
    change_us_K(u1n,u2n,p)
end

function change_us_K(u1,u2,p)
    p.us_K[p.i_state] = u1
    p.us_K[p.i_state+1] = u2
end

# changes the phasis of functions
change_gauge_functions(u1,u2,ξ) = cis(ξ)*u1, cis(-ξ)*u2
    # α = cis(ξ); β = cis(-ξ)
    # px("noRm1 ",norm(u_K(1,p)))
    
    # p.us_K[p.i_state+1] *= β
    # p.us_K[p.i_state+2] *= α
    # px("noRm2 ",norm(u_K(1,p))," α ",α)

# function change_gauge_functions_from_p(ξ,p) # changes the phasis of functions
    # u1, u2 = change_gauge_functions(u_K(1,p),u_K(2,p),ξ,p)
    # change_us_K(u1,u2,p)
    # p.us_K[p.i_state+2] *= α
# end

function test_rot_sym_u(u,c,b;name="u") # tests Ru = cu
    Ru = R_fc_2d(u,b)
    τu = τ_fc_2d(u,shift_K_R(b),b)
    px("Test R "*name*" = c "*name*" : ",distance(Ru,c*τu))
end

function test_rot_sym(p)
    # Tests
    (Ru0,Tu0) = (R_2pi3(u_K(3,p),p),τ_k(u_K(3,p),shift_K_R(p.basis),p)) # u0
    (RS,TS) = (R_2pi3(u_K(1,p),p),τ_k(u_K(1,p),shift_K_R(p.basis),p))   # u1
    (RW,TW) = (R_2pi3(u_K(2,p),p),τ_k(u_K(2,p),shift_K_R(p.basis),p))   # u2

    τau = cis(2π/3)

    px("Test R Φ0 =     Φ0 ",distance(Ru0,Tu0))
    px("Test R Φ1 = τ   Φ1 ",distance(RS,τau*TS)) # R u1 = τ e^{ix⋅(1-R_{2π/3})K} u1 = τ e^{ix⋅[-1,0,0]a^*} u1
    px("Test R Φ2 = τ^2 Φ2 ",distance(RW,τau^2*TW))
    px("Test R v  =      v ",distance(pot_fc(p),R_four(pot_fc(p),p.basis)))
end

function test_M_sym(p)
    (Mu1,Tu1) = (M_change(u_K(1,p),p),τ_k(u_K(1,p),shift_K_M(p.basis),p))   # u1
    (Mu2,Tu2) = (M_change(u_K(2,p),p),τ_k(u_K(2,p),shift_K_M(p.basis),p))   # u2
    (Mu3,Tu3) = (M_change(u_K(3,p),p),τ_k(u_K(3,p),shift_K_M(p.basis),p))   # u3
    px("Test M Φ1 = Φ2 ", distance(Mu1, Tu2))
    px("Test M Φ1 = -Φ2 ",distance(Mu1,-Tu2))
    px("Test M Φ3 = Φ3 ",distance(Mu3,Tu3))
    px("Test M v  =      v ",distance(pot_fc(p),M_four(pot_fc(p),p.basis)))
end

function test_PC_sym(p)
    u1 = u_K(1,p); u2 = u_K(2,p); u3 = u_K(3,p)
    PCu1, PCu2, PCu3 = conj.(u1), conj.(u2), conj.(u3)
    px("Test PC Φ1 = Φ2 ", distance(PCu1, u2))
    px("Test PC Φ3 = Φ3 ", distance(PCu3, u3))
end

function compute_scaprod_Dirac(u1,u2,p)
    vec = ∇_K_common(u2,p)
    c1 = -im*dot(u1,vec[1])
    c2 = -im*dot(u1,vec[2])
    (c1,c2)
end

compute_scaprod_Dirac_u1_u2(p) = compute_scaprod_Dirac(u_K(1,p),u_K(2,p),p)

# Changes the gauge such that <Φ1,(-i∇_x)Φ2> = vF*[1,-i], vF ∈ ℝ+
# the gauge is fixed such that <Φ1,(-i∇_x)Φ2> = <u1,(-i∇_x)u2> = vF*goal where vF ∈ ℝ+, hence <e^(iξ) Φ1,(-i∇_x) e^(-iξ)Φ2> = e^(-2iξ) vF*goal
# this is vF*goal which has to be chosen, and not another one, for our 𝕍 to look like T_BM
function records_fermi_velocity_and_fixes_gauge(u1,u2,p) 
    ψ1, ψ2 = copy(u1), copy(u2)
    gauge_param = 1
    goal = gauge_param*[1;-im]
    (c1,c2) = compute_scaprod_Dirac(ψ1,ψ2,p)
    # px("Values of Fermi velocity params before the gauge change ",c1," ",c2)
    # if !p.alleviate
    # @assert abs(goal[1]/goal[2] - c1/c2) <1e-2
    # end
    I_real = imag(goal[1])<1e-8 ? 1 : 2
    fact = real(goal[I_real])*(I_real==1 ? c1 : c2)
    r,ξ = polar(fact)
    @assert imag(fact*cis(-ξ)) < 1e-10
    ψ1, ψ2 = change_gauge_functions(ψ1, ψ2, ξ/2)

    # Verification that vF is real positive and close to vF
    (c1,c2) = compute_scaprod_Dirac(ψ1,ψ2,p)
    # if !p.alleviate
    x = norm([c1,c2] .- r*goal)
    if x>1e-4
        px("c1,c2=",[c1,c2]," goal=",r*goal," diff=",x)
        px("PROBLEM IN COMMON_GRAPHENE !!! vF IS NOT RIGHT !!! WE DO NOT HAVE GRAPHENE")
    end
    # end
    px("<u1,(-i∇r) u2> = vF[",goal[1],",",goal[2],"] where vF=",r," should be 0.387...")
    # (c11,c22) = compute_scaprod_Dirac(ψ1,ψ2,p)
    # @assert abs(c11-c1)+abs(c22-c2) < 1e-15
    ψ1, ψ2
end

function records_fermi_velocity_and_fixes_gauge_from_p(p) 
    u1 = p.us_K[p.i_state]
    u2 = p.us_K[p.i_state+1]
    ψ1, ψ2 = records_fermi_velocity_and_fixes_gauge(u1,u2,p) 
    change_us_K(ψ1,ψ2,p)
end

function get_u1_u2_in_good_basis(u1,u2,p)
    ψ1, ψ2 = get_natural_basis(u1,u2,p)
    ψ1, ψ2 = records_fermi_velocity_and_fixes_gauge(ψ1,ψ2,p) 
    ψ1, ψ2
end

function get_u1_u2_in_good_basis_from_p(p)
    u1 = p.us_K[p.i_state]
    u2 = p.us_K[p.i_state+1]
    ψ1, ψ2 = get_u1_u2_in_good_basis(u1,u2,p)
    p.us_K[p.i_state] = ψ1
    p.us_K[p.i_state+1] = ψ2
end

function plot_Es_around(p)
    q = 3
    I = max(p.i_state-q,1)
    II = min(p.i_state+q,length(p.Es_K))
    for i=I:II
        eq = i==p.i_state ? "I_STATE=" : ""
        print("i=",eq,i," ",p.Es_K[i],"; ")
    end
end
