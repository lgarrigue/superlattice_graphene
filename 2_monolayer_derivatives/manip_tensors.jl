using LinearAlgebra
px = println
s = sqrt(3)
coef = 96
ω = cis(2π/3)

function is_int(x)
    xx = abs(x)
    y = floor(Int,xx+0.5)
    if abs(xx-y) < 1e-7
        return true
    end
    false
end

# matrix, from numbers to latex
function to_tex(A;name="")
    px("TeXised ",name)
    for i=1:8
        for j=1:8
            a = x2s(A[i,j];tex=true)
            esp = j==8 ? " \\\\" : " & "
            print(a,esp)
        end
        print("\n")
    end
end

function n_ints(A) # number of integers in A
    c = 0
    for i=1:length(A)
        if is_int(A[i])
            c += 1
        end
    end
    c
end

function x2s_real(x;tex=false,error=true,tol=1e-6)
    p = x<0 ? "-" : " "
    xx = abs(x)
    y = floor(Int,xx+0.5)
    ys = floor(Int,xx/s+0.5)
    if abs(xx-y) < tol
        return p*string(y)
    elseif abs(xx-ys*s) < tol
        sa = ys == 1 ? "" : string(ys)
        suf = tex ? "\\sqrt(3) " : "s3"
        return p*sa*suf
    else
        if error
            px("ERREUR IN x2s")
        end
    end
    x
end

function x2s(x;tol=1e-6,tex=false)
    X = real(x)
    Y = imag(x)
    sx = x2s_real(X;tex=tex,tol=tol)
    sy = x2s_real(Y;tex=tex,tol=tol)*"i"
    if abs(X) + abs(Y) < tol
        return "0"
    elseif abs(X) < tol
        return sy
    elseif abs(Y) < tol
        return sx
    end
    sx*"+"*sy
end

get_Q() = [-1 sqrt(3); -sqrt(3) -1]

function get_mats()
    σ1 = [0 1;1 0]
    σ2 = [0 -im;im 0]
    σ3 = [1 0;0 -1]
    D = σ3 + im*σ1
    σ1, σ2, σ3, D
end
# px("Q σ3 Q^T - 4 σ3 D =",Q'*σ3*Q - 4*ω*σ3)
# px("Q D Q^T - 4 y D =",Q'*D*Q - 4*ω*D)

# for α=1:2, β=1:2, γ=1:2
    # px("##### ",α," ",β," ",γ)
    # for a=1:2, b=1:2, c=1:2
        # co = Q[α,a]*Q[β,b]*Q[γ,c]
        # px(a," ",b," ",c," ",x2s(co))
    # end
# end

function get_coords(;dims=3) # dims ∈ [2,3]
    coords_3 = []
    for α=1:2, β=1:2, γ=1:2
        push!(coords_3,(α,β,γ))
    end
    coords_2 = []
    for α=1:2, β=1:2
        push!(coords_2,(α,β))
    end
    dims==2 ? coords_2 : coords_3
end

# coords : [(1, 1, 1), (1, 1, 2), (1, 2, 1), (1, 2, 2), (2, 1, 1), (2, 1, 2), (2, 2, 1), (2, 2, 2)]

function get_W(;dims=3)
    coords = get_coords(;dims=dims)
    Q = get_Q()
    N = 2^dims
    W = [0.0 + 0*im for A=1:N, B=1:N]
    for A=1:N, B=1:N
        cA = coords[A]
        cB = coords[B]
        W[A,B] = prod(Q[cA[j],cB[j]] for j=1:dims)
    end
    # px("Matrix F")
    # display(x2s.(M))
    Matrix(W)
end

P(y) = 4*(y^2 + y + 1)
ze() = zeros(ComplexF64,8,8)+I

function create_C(y)
    C1 = ze(); C2 = ze(); C3 = ze(); C4 = ze(); C5 = ze(); C6 = ze(); C7 = ze(); C8 = ze(); C9 = ze(); C10 = ze(); C11 = ze(); C12 = ze(); C13 = ze(); C14 = ze(); C15 = ze(); C16 = ze(); C17 = ze(); C18 = ze(); C19 = ze(); C20 = ze(); C21 = ze(); C22 = ze(); C23 = ze(); C24 = ze(); C25 = ze(); C26 = ze(); C27 = ze(); C28 = ze(); C29 = ze(); C30 = ze(); C31 = ze(); C32 = ze(); C33 = ze(); C34 = ze(); C35 = ze(); C36 = ze()

    C1[3,:] = [0 -1 1 0 0 0 0 0]
    C2[4,:] = [1 0 0 1 0 0 0 0]
    C3[7,:] = [0 0 0 0 0 -1 1 0]
    C4[8,:] = [0 0 0 0 1 0 0 1]
    C5 = Diagonal([1,1,1/4,1/4,1,1,1/4,-1/4])
    C6[4,:] = [0 0 0 0 0 0 1 0]
    C6[7,:] = [0 0 0 1 0 0 0 0]

    C7[4,:] = [0 0 0 0 1 0 0 0]
    C7[5,:] = [0 0 0 1 0 0 0 0]
    C7[3,:] = [0 0 0 0 0 1 0 0]
    C7[6,:] = [0 0 1 0 0 0 0 0]

    C12[7,:] = [0 0 0 0 0 0 s 2y+1]
    C13[6,:] = [0 0 0 0 2y+1 -s 0 0]
    C14[6,:] = [0 0 0 0 0 0 0 1]
    C14[8,:] = [0 0 0 0 0 1 0 0]

    if abs(y-1)<1e-10
        C16 = Diagonal([1,1,1,1,1,1,1/P(y),1/P(y)])
        C17[5,:] = [0 0 0 0 1 0 0 -3]
        C17[6,:] = [0 0 0 0 0 1 -3 0]
        C18 = Diagonal([1,1,1,1,1/s,1/s,1,1])

        # Bloc du dessus
        C19[1,:] = [1 0 0 0 s 9 -3s -3]
        C19[2,:] = [0 1 0 0 3 s 3 -3s]
        C19[3,:] = [0 0 1 0 3s 3 -s 3]
        C19[4,:] = [0 0 0 1 3 s 3 s]
        C20 = Diagonal([1/2,1/2,1/2,1/2,1,1,1,1])
        C21[2,:] = [s 1 0 0 0 0 0 0]
        C21[3,:] = [-1 0 1 0 0 0 0 0]
        C21[4,:] = [0 1 0 1 0 0 0 0]
        C22[1,:] = [1 -s/4 0 0 0 0 0 0]
        C23 = Diagonal([1/s,1/(4s),1,1,1,1,1,1])

    elseif abs(y-ω)<1e-10
        C16 = Diagonal([1,1,1,1,1/s,1/s,1,1])
        C17[3,:] = [1 0 1 0 0 0 0 0]
        C17[4,:] = [0 -1 0 1 0 0 0 0]
        C18 = Diagonal([1,1,1/(4s),1/(4s),1,1,1,1])

        C19[5,:] = [0 0 0 im 1 0 0 0]
        C19[6,:] = [0 0 -im 0 0 1 0 0]
        C20[6,:] = [0 0 0 0 -im 1 0 0]
        C20[3,:] = [0 0 1 0 1 0 0 0]

        C21[1,:] = [1 0 -im*(3-4s*im) 0 0 0 0 0]
        C21[2,:] = [0 1 s*im 0 0 0 0 0]

        C22[1,:] = [1 0 0 s*im 0 0 0 0]
        C22[2,:] = [0 1 0 im*(3-4s*im) 0 0 0 0]

        C23[1,:] = [1 0 0 0 s 0 0 0]
        C23[2,:] = [0 1 0 0 3 0 0 0]

        C24[2,:] = [-im 1 0 0 0 0 0 0]
        C25 = Diagonal([1,-im/2,1,1,1,1,1,1])
        C26[1,:] = [1 1 0 0 0 0 0 0]
        C27 = Diagonal([1/(-6+2s*im),1/(-s-3*im),1,1,1,1,1,1])

        C28[1,:] = [0 0 1 0 0 0 0 0]
        C28[3,:] = [1 0 0 0 0 0 0 0]
        C28[2,:] = [0 0 0 1 0 0 0 0]
        C28[4,:] = [0 1 0 0 0 0 0 0]

        C29[3,:] = [0 0 0 0 1 0 0 0]
        C29[5,:] = [0 0 1 0 0 0 0 0]

        C30[4,:] = [0 0 0 0 1 0 0 0]
        C30[5,:] = [0 0 0 1 0 0 0 0]

        C31[1,:] = [1 0 0 im 0 0 0 0]
        C31[2,:] = [0 1 0 1 0 0 0 0]
        C31[3,:] = [0 0 1 im 0 0 0 0]
    else
        px("y=",y," PAS PRIS EN COMPTE, DONC PAS CENSE MARCHER !!!!!!!!!!!!!!!!!!!!!!")
    end

    C = C35*C34*C33*C32*C31*C30*C29*C28*C27*C26*C25*C24*C23*C22*C21*C20*C19*C18*C17*C16*C15*C14*C13*C12*C11*C10*C9*C8*C7*C6*C5*C4*C3*C2*C1
    coef*C
end

function get_tilde(y;verbose=false)
    # COMPUTING TILDE
    Cstr = "B"
    if abs(y-1) < 1e-10
        Cstr = "A"
    end
    print("\nCOMPUTE ",Cstr," TILDE ")

    # Compute W
    W = get_W(;dims=3)
    Wm8 = W-8*y*I
    # Compute A or B
    C = create_C(y) # C = A or B
    if y==1
        @assert norm(imag.(C))<1e-13
        C = real.(C)
    end
    # invC = inv(C)
    # px("C^(-1) C ")
    # display(x2s.(invC*C))
    # px("C^(-1) ")
    for j=1:8
        # px("invC[:,",j,"]")
        # display(x2s.(invC[:,j]))
        # display(x2s.(invC[:,j]))
    end
    # px("number of integers ",100*n_ints(C)/length(C),"%")
    if verbose
        px("\n \n ",Cstr," =")
        display(x2s.(C))
        # display(C)
    end
    px("y=",y," ; det(",Cstr,") =",det(C))
    C_tilde = C*Wm8/coef
    C, C_tilde
end

function nzeros(A)
    nz = 0
    for m in A
        if abs(m)<1e-3
            nz+=1
        end
    end
    nz
end

function do_stuff()
    W = get_W(;dims=3)
    A, Atilde = get_tilde(1;verbose=true)
    B, Btilde = get_tilde(ω;verbose=true)

    # for N in [Atilde,Btilde]
        # nze = nzeros(N)
        # px("Number of zeros ",nze," ",100*nze/8^2,"%")
    # end
    to_tex(get_W(;dims=3);name="W")
    to_tex(Atilde;name="Atilde")
    to_tex(Btilde;name="Btilde")
end

# f.(M)
do_stuff()
