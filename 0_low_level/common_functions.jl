# using Plots#, LaTeXStrings
using DelimitedFiles, FFTW
px = println

import Base.+  
+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)  

new_fun_2d(b) = zeros(ComplexF64,b.N,b.N)
new_fun_3d(b) = zeros(ComplexF64,b.N,b.N,b.Nz)

####################################### Miscellaneous low level stuff

hartree_to_ev = 27.2114
ev_to_hartree = 1/hartree_to_ev
fill1d(x,n) = [deepcopy(x) for i=1:n]
fill2d(x,n) = [deepcopy(x) for i=1:n, j=1:n]
rotM(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]
polar(x) = abs(x),atan(imag(x),real(x))
cyclic_conv(a,b,Vol) = myfft(myifft(a,Vol).*myifft(b,Vol),Vol)/sqrt(Vol) # convolution of Fourier vectors. DON'T USE IT FOR MULTIPLICATIONS WITH THE /sqrt(Vol)

function distance(f,g) # relative distance between f and g, whatever they are
    nor = 0.5*(norm(f) + norm(g))
    if abs(nor) < 1e-15
        px("Division by zero in distance")
        return 0
    end
    norm(f.-g)/nor
end

# if the directory path doesn't exist, it creates it, otherwise does nothing
create_dir(path) = if !isdir(path) mkdir(path) end

function re(x)
    @assert imag(x) < 1e-10
    x = real(x)
end
