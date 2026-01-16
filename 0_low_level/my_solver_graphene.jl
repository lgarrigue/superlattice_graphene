# Schrodinger and its resolution

function diag_monolayer_at_k_mine(kred,v_dir,n_eigenvectors,b;tol=1e-4,maxiter=100,full_diag=false,coef_kin=1,prints=false)
    L = action_for_k(kred,v_dir,b;coef_kin=coef_kin)
    kns = X("kcarts_sq",b.K_red_2d,b.stock)
    @assert length(v_dir) == b.Nfull
    precon = vcat(kns...)
    n_max = b.Nfull
    full_diag2 = full_diag || n_max < 100
    (λs,cv,Xf) = solve_lobpcg(L,n_eigenvectors,precon; maxiter=maxiter,tol=tol,X0=-1,full_diag=full_diag2,prints=prints)
    @assert size(Xf,2) == n_eigenvectors
    # Extracts eigenvectors
    ref_vec = ones(b.Nfull)
    ref_gauge = ref_vec/norm(ref_vec) 
    n_eigenmodes = length(λs)
    us = []
    for i=1:n_eigenvectors
        u_lin = Xf[:,i]
        # Fixes the gauges
        λ = sum(conj(u_lin).*ref_gauge)
        c = λ/abs(λ)
        u_lin *= c
        push!(us,u_lin)
    end
    (λs, us)
end
