include("../1_common_graphene_manips/common_graphene.jl")

function verifications()
    # N = 27; Nz = 600; gauge_param = 1
    # N = 27; Nz = 160; gauge_param = 1
    # N = 20; Nz = 120; gauge_param = 1
    p = get_standard_param_module("dftk_import")

    px("Generates Hamiltonian")
    get_dirac_eigenmodes(p;check_K_cart=false)
    px("Finished generating Hamiltonian")

    # Verifications. Don't do them when the size of matrices is large, this will crash
    ξ = ones(ComplexF64,N_basis(p))
    # Checks that applying p.spec.H_K_block*ξ and p.spec.H_matrix are equivalent
    @assert distance(p.spec.H_K_block*ξ, p.spec.H_matrix*ξ) < 1e-13

    # Verifications on operators
    px("Number of blocks ",length(p.spec.ham.blocks))
    ops = p.spec.ham.blocks[1].operators
    fourier_ops  = filter(o -> o isa DFTK.FourierMultiplication,   ops)
    real_ops     = filter(o -> o isa DFTK.RealSpaceMultiplication, ops)
    nonlocal_ops = filter(o -> o isa DFTK.NonlocalOperator,        ops)
    divAgrad_ops = filter(o -> o isa DFTK.DivAgradOperator,        ops)
    strs = ["total","fourier","real","nonloc","divAgrad"]
    terms = [ops,fourier_ops,real_ops,nonlocal_ops,divAgrad_ops]
    px("Number of operators ")
    for i=1:length(strs)
        n = length(terms[i])
        print(" ",strs[i]," ",n,"\n")
    end


    # Non local
    E_non_local = 0
    nlV = -1
    for i=1:length(nonlocal_ops)
        nlV = Matrix(nonlocal_ops[i])
        E_non_local += re(dot(u_K(1,p),nlV*u_K(1,p)))
    end
    @assert E_non_local < 1e-12
    if E_non_local > 1e-9
        px("Energy nonlocal ",E_non_local)
    end

    # Potentials
    E_pot_tot = 0
    v_summed = 0*p.v_dir
    u1_dir = u_dir_n(1,p)
    Vtot = 0*nlV
    for i=1:length(real_ops)
        op = real_ops[i]
        V = Matrix(op)
        Vtot += V
        E = re(dot(u_K(1,p),V*u_K(1,p)))
        E_pot_tot += E
        # px("Energy potential from Fourier ball ",i," ",E)
        v = op.potential
        Edir = re(scaprod(u1_dir,v.*u1_dir,p.basis,false))
        @assert distance(E,Edir) < 1e-10
        v_summed .+= v
    end
    
    v_tot = DFTK.total_local_potential(p.spec.ham)[:,:,:,1]
    @assert distance(v_tot,p.v_dir) < 1e-10
    @assert distance(v_summed,p.v_dir) < 1e-10
    # px("Total potential energy ",E_pot_tot)

    # Potential energy
    uv = mul_3d(pot_fc(p),u_fcl_n(1,p),p.basis)
    # uv = mul_3d(u_fcl_n(1,p),K_red_2d pot_fc(p), p.basis)
    pot_en_C = scaprod(u_fcl_n(1,p),uv,p.basis)
    @assert imag(pot_en_C) < 1e-12
    pot_en = re(pot_en_C)
    pot_dir = re(scaprod(u1_dir, p.v_dir.*u1_dir, p.basis,false))
    @assert distance(pot_dir,pot_en) < 1e-10
    @assert distance(E_pot_tot,pot_en) < 1e-10 # comparision with the fb one

    # Potential matrix
    ψ = Vtot*[i==1 ? 1 : 0 for i=1:N_basis(p)]
    x = distance(Vtot[:,1],ψ)
    px("Test on matrix potential ",x)

    # Verifications on Fermi energy
    Ef1 = re(dot(u_K(1,p), p.spec.H_K_block * u_K(1,p)))
    Ef2 = re(dot(u_K(2,p), p.spec.H_K_block * u_K(2,p)))

    E1s = re(p.Es_K[p.i_state])
    E2s = re(p.Es_K[p.i_state+1])
    @assert sum([abs(Ef1-x) for x in [Ef2,E_fermi(p),E1s,E2s]]) < 1e-5
    px("Fermi energy ",Ef1)

    # Verification on kinetic
    kin_op1 = fourier_ops[1]
    # kin_op2 = DFTK.Kinetic()
    # Kψ2 = kin_op1 .* u_K(1,p)
    Kψ1 = kin_op1.multiplier .* u_K(1,p)
    kin_en_mul = re(dot(u_K(1,p),Kψ1))
    # px("dist ",norm(Hψ1 .- Hψ2)/norm(Hψ1))
    # kin_en1 = scaprod(u,Hb*u,p.basis)
    kin_en = -0.5*re(scaprod(u_fcl_n(1,p), Δ_k_mine(u_fcl_n(1,p),p.basis;kred=p.basis.K_red_3d),p.basis))
    kin_en_fb = -0.5*re(dot(u_K(1,p), Δ_K_dftk(u_K(1,p),p)))
    px("kin_en_mul ",kin_en_mul," kin_en ",kin_en," kin_en_fb ",kin_en_fb)
    @assert distance(kin_en, kin_en_mul) < 1e-9
    @assert distance(kin_en, kin_en_fb) < 1e-9

    # Total energy
    Etot_schro = u_fcl_n(1,p)⋅apply_schro(u_fcl_n(1,p),p.v_dir,p.basis;kred=p.basis.K_red_2d)
    @assert imag(Etot_schro) < 1e-11
    Etot_schro = real(Etot_schro)
    @assert distance(Etot_schro,E_fermi(p)) < 1e-9

    Etot = kin_en+pot_en
    @assert distance(Etot,E_fermi(p)) < 1e-9
    px("Kinetic energy: ",kin_en) # 1.30156
    px("Potential energy: ",pot_en) # -1.3620
    px("Total energy: ",Etot) # -0.06046 = E_fermi

    # Verifications on residual
    resi_fb = p.spec.H_K_block * u_K(1,p) -  Ef1 * u_K(1,p)
    res = norm(resi_fb)
    # @assert res < 1e-4
    px("Verification residual Fourier ball ",res)
    mΔu = -0.5*Δ_k_mine(u_fcl_n(1,p),p.basis;kred=p.basis.K_red_3d)
    resi_fc = mΔu + uv - E_fermi(p)*u_fcl_n(1,p)
    res_fc = norm(resi_fc)
    px("Verification residual Fourier cube ",res_fc)
    # It is bad, but after a study, only the relatively high modes of res_fc are high, so res_fc is small on the fourier modes close to zero
    c = 0
    for u in p.us_K
        ufc = u_fcl(u,p)
        c += abs(ufc⋅resi_fc)
    end
    px("Scalare residual with other vectors ",c) # is small
    hu_fc1 = mΔu + uv
    hu_fc2 = u_fcl(p.spec.H_K_block * u_K(1,p), p)
    px("Norms ",norm(hu_fc1)," ",norm(hu_fc2))
    hu_fc1 /= norm(hu_fc1)
    hu_fc2 /= norm(hu_fc2)
    px("Compare commutation, equality ?",norm(hu_fc1-hu_fc2)) # doesn't work

    # Verification on pseudo-inverse
    px("(H-E) u_K_fb = 0 ? ",res)
    u = ones(ComplexF64,N_basis(p))
    Ru = apply_pseudo_inverse(u,p)
    px("norm R u=",norm(Ru),", u being only ones")
    di = distance(p.Es_K[p.i_state]*Ru - p.spec.H_K_block*Ru, u)
    px("(H-E)R=1 ? ",di)
    di = norm(apply_pseudo_inverse(p.us_K[p.i_state],p))
    px("R u_K = 0 ? ",di)

    # pot_real = DFTK.irfft(p.basis, pot_fc(p))
    # pot_op = DFTK.RealSpaceMultiplication(p.basis, p.K_kpt, pot_fc(p))

    # mat_kin = Matrix(kin_op)
    # mat_pot = Matrix(pot_op)

    # operators = [kin_op]
    # Hb = HamiltonianBlock(p.basis, p.K_kpt, operators; scratch=nothing)
    nothing
end

verifications()
