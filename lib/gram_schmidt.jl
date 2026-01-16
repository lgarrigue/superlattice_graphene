using Random # for shuffle

function gram_schmidt_one_vector(v,V)
    if length(V)==0
        return v/norm(v)
    end
    u = v/norm(v)
    u_perp = u
    u_perp -= sum((vk⋅u)*vk for vk in V) 
    # px("norm of added vector ",norm(u_perp))
    u_perp /= norm(u_perp)
    u_perp
end

function gram_schmidt_family_without_forcing(V)
    new_V = [V[1]/norm(V[1])]
    n = length(V)
    for j=2:n
        u_perp = gram_schmidt_one_vector(V[j],new_V[1:j-1])
        push!(new_V,u_perp)
        # test_orthonormal_basis(new_V;tol = 1e-8)
        # px("TEST ORTH ",j)
    end
    # test_orthonormal_basis(new_V;tol = 1e-8)
    new_V
end

function gram_schmidt_family(V;force=false,tol=1e-10,print=false)
    l = gram_schmidt_family_without_forcing(V)
    err = orthonormality(l)
    if !force
        if print px("Orth error after Gram-Schmidt : ",err) end
        return l
    end
    max_iter = 50
    j = 0
    while err > tol && j < max_iter
        l = gram_schmidt_randomly(l;n=1)
        err = orthonormality(l)
        j += 1
    end
    if print || err > tol
        px("Gram-Schmited ",j," times more")
        px("Ortho error after Gram-Schmidt=",err)
    end
    if err > tol
        px("GRAM SCHMIDT NOT SUCCESSFULL after max number of iterations, error is ",err)
    end
    l
end

# dont_touch : indices of the ones we do not modify
function gram_schmidt_randomly(l;n=1,dont_touch=[]) 
    len = length(l)
    ll = copy(l)
    for j=1:n
        nl = []
        s = shuffle((1:len))
        l_interm = [ll[k] for k in s]
        for i=1:len
            u_perp = gram_schmidt_one_vector(l_interm[i],nl)
            push!(nl,u_perp)
        end
        ll = copy(nl)
    end
    ll
end

function orthonormality(l)
    N = length(l)
    S = [l[i]⋅l[j] for i=1:N, j=1:N]
    norm(S-I)/norm(S)
end

function test_orthonormal_basis(l;tol = 1e-8)
    N = length(l)
    err = orthonormality(l)
    px("test orthonormal basis ",err," N=",N)
    @assert err < tol
end
