# were previously in commun_functions.jl

####################################### Operations on functions in Fourier space
function myfloat2int(x;warning=true,name="") # float to int
    y = floor(Int,x+0.5)
    if abs(x-y)>1e-5 && warning
        px("NOT AN INT in ",name," ",x," ",y)
    end
    y
end

# translation u(x) := u(x - ya) = ∑_m u_m e^{ima^*(x-ya)}/sqrt(Vol) has Fourier coefs u_m e^{-i2π m⋅y}
function translation_interpolation(u_f,y_red,b) # y_red in reduced direct space, u_f in Fourier
    dim = length(size(u_f))
    if dim==3
        return [u_f[i,j,l]*cis(-2π*y_red⋅[b.kred_x[i],b.kred_x[j]]) for i=1:b.N, j=1:b.N, l=1:b.Nz]
    else
        return [u_f[i,j]*cis(-2π*y_red⋅[b.kred_x[i],b.kred_x[j]]) for i=1:b.N, j=1:b.N]
    end
end

function R_four(a,b) 
    R_four_2d = matrix_rot_red(-2π/3,b)
    apply_map_four(X ->  R_four_2d*X,a,b) # rotation of 2π/3, in Fourier space
end

function M_four(a,b)
    M_four_2d = float2int.(cart2red_mat([1 0;0 -1],b);name="M") # mirror symmetry, M_four_2d = [0 1;1 0]
    apply_map_four(X ->  M_four_2d*X,a,b) # Mu(x,y) := u(x,-y), mirror
end

# J_four(a,b) = apply_map_four(X -> -[1 -2;2 -1]*X,a,b) # rotation of -π/2, in Fourier space, with scaling of sqrt(3)
parity_four(a,b) = apply_map_four(X -> -X,a,b)
conj_four(a,b) = apply_map_four(X -> -X,conj.(a),b) # if f is in direct space, and g(x) := conj.(f(x)), then hat(g)_m = conj(hat(f))_{-m}
σ1_four(a,b) = apply_map_four(X -> [0 1;1 0]*X,a,b)

function apply_map_four(L,u,b) # gives C_m = u_{Lm}, ie applies L in Fourier reduced space. u is in Fourier
    dim = length(size(u)); tu = typeof(u[1])
    a = dim==2 ? zeros(tu,b.N,b.N) : zeros(tu,b.N,b.N,b.Nz)
    for K=1:b.N, P=1:b.N
        k0 = b.kred_x[K]; p0 = b.kred_x[P]
        c = L([k0,p0])
        k1 = c[1]; p1 = c[2]
        (k2,p2) = k_inv(k1,p1,b)
        if dim==2
            a[K,P] = u[k2,p2]
        else 
            a[K,P,:] = u[k2,p2,:]
        end
    end
    a
end

function apply_map_four_back(L,u,b) # gives C_m such that C_{Lm} = u_m
    a = zeros(typeof(u[1]),b.N,b.N)
    for K=1:b.N, P=1:b.N
        k0 = b.kred_x[K]
        p0 = b.kred_x[P]
        c = L([k0,p0]); k1 = c[1]; p1 = c[2]
        (k2,p2) = k_inv(k1,p1,b)
        a[k2,p2] = u[K,P]
    end
    a
end

function test_nΩ_periodicity(u,n,b) # test whether ∀j u(x + 1/n ja) = u(x) ⟺ ∀j,m u_m(e^{ijm/n}-1)=0 ⟺ ∀j,m jm ∈nℕ+1 ⟹ u_m = 0 ⟺ ∀m m1 mod n ≠0 or m2 mod n ≠0 ⟹ u_m=0
    c = 0
    for ik_m=1:b.Nfull
        k = ik2kred_2d(ik_m,b)
        if mod(k[1],n)!=0 || mod(k[2],n)!=0
            c += abs(u[ik_m])
        end
    end
    c /= norm(u)
    px("Is ",n,"Ω-periodic ? ",c)
end

J_four_back(a,b) = apply_map_four_back(X -> -[1 -2;2 -1]*X,a,b) # rotation of -π/2, in Fourier space, with scaling of sqrt(3)

####################################### Coordinates changes in periodic direct space
# Example
# xi = [1 2 3 4]
#  x = [0 1 2 3]
xi2x(x,b) = mod(x-1,b.N)
x2xi(x,b) = mod1(x+1,b.N)

# Xi = xi + n*N in natural coordinates corresponding to xi
# Xpi = (xi,n)
# X = Xi-1
# Xred = x_red[xi] + n*b.L
# Xp = (x,n)
Xpi2Xp(Xpi,b) = (xi2x(Xpi[1],b),Xpi[2])
Xp2Xpi(X,b) = (x2xi(X[1],b),n)
Xpi2Xred(X,b) = b.x_axis_red[X[1]] + b.La0*X[2] # NECESSITATES THAT L1 == L2 !!!

Xi2Xpi(xi,b) = (mod1(xi,b.N),fld1(xi,b.N))
X2Xi(X) = X+1
Xi2X(Xi) = Xi-1
X2Xpi(X,b) = Xi2Xpi(X2Xi(X),b)

####################################### Operations on functions in direct space

# Applies coordinates transformation M
function apply_coordinates_operation_direct(M,b) # (Op f)(x) = f(Mx), 2d or 3d
    function f(ϕ,b)
        ψ = similar(ϕ)
        dim = length(size(ϕ))
        for xi=1:b.N, yi=1:b.N
            X = Int.(M([xi2x(xi,b);xi2x(yi,b)]))
            Xi = x2xi(X[1],b); Yi = x2xi(X[2],b)
            if dim == 3
                ψ[xi,yi,:] = ϕ[Xi,Yi,:]
            else
                ψ[xi,yi] = ϕ[Xi,Yi]
            end
        end
        ψ
    end
    f
end

function apply_coordinates_operation_direct_on_z(M,b)
    function f(ϕ,b)
        ψ = similar(ϕ)
        dim = length(size(ϕ))
        for zi=1:b.N
            z = Int.(M(xi2x(zi,b)))
            zi = x2xi(z,b)
            ψ[:,:,zi] = ϕ[:,:,zi]
        end
        ψ
    end
    f
end

# R(ϕ,b)             = apply_coordinates_operation_direct(X -> R_four_2d*X,b)(ϕ,b)
translation(ϕ,v,b) = apply_coordinates_operation_direct(X -> X.-v,b)(ϕ,b)
parity_x(ϕ,b)      = apply_coordinates_operation_direct(X -> -X,b)(ϕ,b)

translation2d(u,a,b) = [u[mod1(x-a[1],b.N),mod1(y-a[2],b.N)] for x=1:b.N, y=1:b.N]

function parity_z(u,b)
    u_fc = myreshape(u,b)
    u_p = ""
    if b.dim == 3
        u_p = [u_fc[x,y,mod1(2-z,b.Nz)] for x=1:b.N, y=1:b.N, z=1:b.Nz]
    else
        u_p = [u_fc[mod1(2-z,b.Nz)] for z=1:b.Nz]
    end
    mylinearize(u_p)
end

# Applies coordinates transformation M and does a Bloch transform at the same time
# (Op B f)(x) = (Op ∘ exb)(x) * (Op ∘ f)(x) = exp(Mx) * f(Mx)
function apply_Op_B(M,k,b) 
    function f(u,k,b)
        u_fc = myreshape(u,b)
        ψ = similar(u_fc)
        for xi=1:b.N, yi=1:b.N
            # Rotation
            RX = Int.(M([xi2x(xi,b);xi2x(yi,b)]))
            Xpi = X2Xpi(RX[1],b); Ypi = X2Xpi(RX[2],b)

            # Phasis
            X = Xpi2Xred(Xpi,b); Y = Xpi2Xred(Ypi,b)
            φ = cis(2π*k⋅[X;Y])
            ψ[xi,yi,:] = φ*u_fc[Xpi[1],Ypi[1],:]
        end
        mylinearize(ψ)
    end
    f
end

function OpL_fc_2d(u, L, b) 
    Lu = zero(u)
    for ik=1:b.Nfull
        k = ik2kred_2d(ik,b)
        Lk = L(k)
        ik_L = kred2ik_2d(Lk,b)
        if ik_L==-1
            Lu[ik] = 0
        else
            Lu[ik] = u[ik_L]
        end
    end
    Lu
end

function R_fc_2d(u,b)
    R_four_2d = matrix_rot_red(-2π/3,b)
    L(k) = R_four_2d*k
    OpL_fc_2d(u,L,b)
end

function M_fc_2d(u,b)
    M_four_2d = float2int.(cart2red_mat([1 0;0 -1],b);name="M") # mirror symmetry, M_four_2d = [0 1;1 0]
    L(k) = M_four_2d*k
    OpL_fc_2d(u,L,b)
end

function P_fc_2d(u,b)
    P_four_2d = float2int.(cart2red_mat([-1 0;0 -1],b);name="P") # parity
    L(k) = P_four_2d*k
    OpL_fc_2d(u,L,b)
end

function τ_fc_2d(u,kred,b)
    OpL_fc_2d(u,k -> k .- kred[1:2],b)
end

function test_R_inv(u_fc,b;name="") # ! in Fourier, not in direct space
    R_u_fc = R_fc_2d(u_fc,b)
    x = norm(R_u_fc - u_fc)/norm(u_fc)
    px("R "*name*" = "*name*" ",x)
    x
end

# Bloch transform and Rotation, RBu = Rexp(...) * Ru. We obtain (RBu)(x,y) for (x,y) ∈ x_grid_red
function RB(u,k,b)
    R_four_2d = matrix_rot_red(-2π/3,b)
    apply_Op_B(X -> R_four_2d*X,k,b)(u,k,b)
end

PB(u,k,b) = apply_Op_B(X -> -X,k,b)(u,k,b)

####################################### Test functions

function test_scaprod_fft_commutation(b) # tests commutation between FFT and scalar product operations
    ϕ = randn(b.N,b.N,b.N)
    ψ = randn(b.N,b.N,b.N)
    c = scaprod(ϕ,ψ,p,false) - scaprod(myfft(ϕ,b.Vol),myfft(ψ,b.Vol),b)
    px("Test scalar product ",c)
end

function test_k_inv() # tests kz inversion
    N = 10
    kred_x = fftfreq(N)*N
    c = true
    for ki=1:N
        k = kred_x[ki]
        if kred_x[k_inv_1d(k,N)] != k
            c = false
        end
    end
    px("Test k inv: ",c ? "good" : "problem")
end

test_hermitianity(M,name="") = px(string("Test Hermitianity ",name," : "),antiherm(M)) # tests M^*=M

function test_x_parity(u,p;name="") # Tests u(-x) = u(x) (or u(-x,z) = u(x,z)), where u is in direct space
    c = sum(abs.(parity_x(u,b) .- u))/sum(abs.(u))
    px("Test ",name,"(-x) = ",name,"(x) : ",c)
end

# tests u(-z) = ε u(z)
test_z_parity(u,ε,b;name="function") = px("Test ",name,"(-z) = ",ε==-1 ? "-" : "",name,"(z) : ",norm2(u.- ε*parity_z(u,b.basis),b)/norm2(u,b.basis))

intZ(f,b) = b.dz*[sum(f[x,y,:]) for x=1:size(f,1), y=1:size(f,2)] # partial integration over XY
intXY(f,b) = b.dS*[sum(f[:,:,z]) for z=1:size(f,3)] # partial integration over Z
average_over_xy(f_dir,b) = [sum(f_dir[:,:,z]) for z=1:b.Nz]/b.N^2 # partial integration over Z

z_translation(a,Z,b) = [a[x,y,mod1(z-Z,b.Nz)] for x=1:b.N, y=1:b.N, z=1:b.Nz] # translation on Z
r_translation(a,s,b) = [a[mod1(x-s[1],b.N),mod1(y-s[2],b.N),z] for x=1:b.N, y=1:b.N, z=1:b.Nz] # translation on XY, s ∈ {0,…,b.N-1}^2, 0 for no translation

matrix_rot_red(θ,b) = float2int.(cart2red_mat(rotM(θ),b);name="mat_rot") # rotation operator in reduced Fourier space

function test_rotations(period;first=[1,0])
    n = period
    N = 13
    b = graphene_basis_2d(N;a=3)
    R = matrix_rot_red(2π/n,b)
    set = [(R^i)*first for i=0:n]
    px("Set of R 2π/",n," : ",set)
end
