# include("../0_low_level/operations_on_fc_functions.jl")
function generate_some_sym_potential_dir(V_intensity,b;fourier=false,first_k=[1,0],mykreds=[])
    kreds = mykreds
    if kreds==[]
        R_four_2d = matrix_rot_red(-2π/3,b)
        kreds = [(R_four_2d^i)*first_k for i=0:2] # corresponds to Fefferman-Weinstein condition. Do not put [first_k,first_k], we are not in the same basis as them
        kreds = vcat(kreds,-kreds) # need to have e^{imax} and e^{-imax} to have a real potential, so there are only cosinuses
    end

    # px("ROT ",[(matrix_rot_red(2π/3,b)^i)*first_k for i=0:3])
    px("KREDS ",kreds)
    # px("R_four = ",R_four_2d)
    # kreds = [[1,0],[1,1]]
    # px("kreds ",kreds)
    v_fc = new_fun_2d(b)
    for kred in kreds
        ik = kred2ik_2d(kred,b)
        if 1 ≤ ik ≤ length(v_fc)
            v_fc[ik] = 1 # when it's positive, band diagram looks like graphene bands more, don't know why, but it's equivalent a priori
        else
            px("Problem producing potential_dir in my_structs_2d, increase size of basis")
            @assert false
            return -1
        end
    end
    v_fc *= V_intensity
    if fourier
        return v_fc
    end
    myifft(v_fc,b.Vol)
end

function generate_some_2d_param(N,V_intensity;n_eigenmodes=-1,first_k=[1,0],i_state=1,a=5)
    # Basis
    b = Basis(2)
    b.N = N
    init_cell_vectors(b;a=a) 
    init_cell_infinitesimals(b)

    # Parameters
    p = Params(2)
    p.basis = b

    # Potential
    p.v_dir = generate_some_sym_potential_dir(V_intensity,b;first_k=first_k)

    # Fills eigenmodes at K
    get_dirac_eigenmodes_mine(p;n_eigenmodes=n_eigenmodes)

    #Other
    p.i_state = i_state
    p.type = "mine_2d"

    # Rotate eigenmodes in the right way
    get_u1_u2_in_good_basis_from_p(p)
    # test_R_inv(myfft(p.v_dir,p.basis.N),p.basis;name="v_dir")
    test_rot_sym(p)
    
    p
end

function get_dirac_eigenmodes_mine(p;n_eigenmodes=-1)
    # n_eigenmodes = p.basis.Nfull
    # Do the diagonalization
    l = p.basis.Nfull
    if n_eigenmodes!=-1 && n_eigenmodes<l
        l = n_eigenmodes
    end
    K = p.basis.K_red_3d
    (p.Es_K, p.us_K) = diag_monolayer_at_k_mine(K,p.v_dir,l,p.basis;tol=1e-6)
    (p.Es_K, p.us_K)
end
