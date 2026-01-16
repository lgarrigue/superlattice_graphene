# Tests if a family of vector is linearly independent or not
using LinearAlgebra#, MKL

# l is a list of vectors or same size
function test_linear_independence(l;tol=1e-3)
    n = length(l)
    if n > length(l[1])
        px("PROBLEM IN test_linear_independence, lengths")
    end
    A = [l[j][i] for i=1:length(l[1]), j=1:length(l)]
    Q, R = qr(A)
    px("Sizes R ",size(R,1)," ",size(R,2))
    @assert n == size(R,1)
    @assert n == size(R,2)
    diag = [R[i,i] for i=1:n]
    free = true
    print("Up to tolerence ",tol," : ")
    for i=1:n
        x = abs(diag[i])
        if x<tol
            print("l is not free, i=",i," qr=",x," ")
            free = false
        end
    end
    if free
        print("l is free, min QR diag is ",minimum(abs.(diag)))
    end
    # print("\nFull diag is ",diag)
    print("\n")
end

function test_of_test_linear_independence()
    l1 = [1,2,-1,7,4,6]
    l2 = [1,-1,-3,4,1,2]
    l3 = [0,1,-7,1,1,2]
    l4 = l1 + 2*l2 + 3*l3
    l_not_free = [l1,l2,l3,l4]
    l_free = [l1,l2,l3]

    test_linear_independence(l_not_free;tol=1e-3)
    test_linear_independence(l_free;tol=1e-3)
end

# test_of_test_linear_independence()
