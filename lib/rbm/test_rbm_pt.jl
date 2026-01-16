include("../pt/derivatives_stock.jl")
include("../pt/asymptotic_behavior.jl")
include("rbm_module.jl")
# include("../../rbm_pt/schrodinger.jl")

herm(H) = (H+H')/2
random_Hn(N) = [copy(herm(randn(N,N))) for k=0:1]
dist(a,b) = norm(a.-b)/norm(a)

function schro_Hn(N)
    (A,H0,G0,G1,G2) = create_schrodinger_potentials(N;plot=false)
    [H0 .+ G0,G1]
end

function compute_at_H(H,R)
    E_exact, ϕ_exact = eigenmodes(H)
    E_rbm, ϕ_rbm = rbm_approx_eigenmodes(H,R)
    ψ = ϕ_exact[1]
    φ = ϕ_rbm[1]
    dist_quo(ψ,φ)
end

function compute_results(Hλ,R)
    test_orthonormal_basis(R)
    λs = get_λs(1,-1.2) # prepares some λ
    ds = []
    # Computes solutions
    for λ in λs
        di = compute_at_H(Hλ(λ),R)
        push!(ds,di)
    end
    ds, λs
end

function test()
    N = 30; ℓ = 3
    # Hs = random_Hn(N)
    Hs = schro_Hn(N)
    S = stock_at_H(Hs) # computes derivatives
    R = EigenRBM(ℓ,0) # builds RBM basis
    init_R_from_first_derivatives(R,S)
    do_Gram_Schmidt(R;n=20)
    test_R(S,R)
    Hλ(λ) = Hs[1] + λ*Hs[2]
    ds, λs = compute_results(Hλ,R)
    print_results(λs,ds,ℓ+1)
end

# test()
