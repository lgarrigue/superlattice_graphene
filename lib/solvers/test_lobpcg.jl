include("../../lib/solvers/lobpcg.jl")

function test_lobpcg()
    N = 1000
    A = randn(ComplexF64,N,N)
    A = Hermitian((A+A')/2)
    L = X -> A*X
    precon = ones(N)
    # precon = zeros(N) # with this one, does not work !
    full_diag = false
    tol = 1e-4
    n_eigenvectors = floor(Int,N/10)
    println("Starts LOBPCG")
    (λs,cv,Xf) = solve_lobpcg(L,n_eigenvectors,precon; maxiter=100,tol=tol,X0=-1,full_diag=full_diag)
    println("LOBPCG Converged ? ",cv)
end

test_lobpcg()
