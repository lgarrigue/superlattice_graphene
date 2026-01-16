######################### Operations on Fourier ball functions for DFTK manipulation



# from a list G=[a,b,c] of int, gives the iG such that Gvectors[iG]==G. If G is not in it, gives nothing
function index_of_Gvector(G,p) 
    tupleG = Tuple(G)
    if haskey(p.spec.Gvectors_inv,tupleG)
        return p.spec.Gvectors_inv[tupleG]
    else
        return nothing
    end
end

# Generates the operator doing U_m -> V_m := U_{Lm}
# L is an action on the reduced Fourier space
function OpL(u, p, L) 
    Lu = zero(u)
    for iG=1:N_basis(p)
        FiG = p.spec.Gvectors[iG]
        LFiG = L(FiG)
        FinvLFiG = index_of_Gvector(LFiG,p)
        if isnothing(FinvLFiG)
            Lu[iG] = 0
        else
            Lu[iG] = u[FinvLFiG]
        end
    end
    Lu
end

M2d_2_M3d(M) = [M           zeros(2,1);
                zeros(1,2)          1 ]

# (Ru)_G = u_{M G} or (Ru)^D_m = u^D_{F^{-1} M F(m)}, where ^D means that we take the discretization vector
function R_fb(u,p)
    R_four_2d = matrix_rot_red(-2π/3,p.basis)
    L(G) = M2d_2_M3d(R_four_2d)*G
    OpL(u,p,L)
end

function M_fb(u,p) # mirror
    M_four_2d = float2int.(cart2red_mat([1 0;0 -1],p.basis);name="M")
    L(h) = M2d_2_M3d(M_four_2d)*h
    OpL(u,p,L)
end

# τ is equivalent to multiplication by e^{i cart(k) x), where k is in reduced, e^{i m0 a^*⋅x} u = ∑ e^{ima^*⋅x} u_{m-m0}
τ_fb(u,k,p) = OpL(u,p,G -> G .- k)
P_fb(u,p) = OpL(u,p,G -> -G) # parity

inverse_dict_from_array(a) = Dict(value => key for (key, value) in Dict(zip((1:length(a)), a)))
