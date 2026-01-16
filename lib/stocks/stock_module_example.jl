using CairoMakie 
include("stock_module.jl")
px = println

# all the parameters of one computation
mutable struct OneComputation
    n
    k
    color
    function OneComputation()
        oc = new()
        oc.n = 0
        oc.k = 0
        oc.color = "red"
        oc
    end
end

function create_stock()
    # 1) Define the different kinds of recursive objects we will need, just for clarity here
    # kinds = ["z","β"]

    # 2) Define default parameters values
    p = Dict()
    p["y"] = 1

    # 3) Define how to compute the quantities
    function compute_z(coords,St::Stock)
        n = coords
        x = -1
        if n == 1
            x = 1
        else
            x = 1/factorial(big(n)) + St.p["y"]*sum([X("z",s,St)/factorial(big(n-s)) for s=1:n-1])
        end
        px("Computed z(n)=z(",n,") = ",x)
        x
    end

    function compute_β(coords,St::Stock)
        (n,k) = coords
        x = 0
        if k >= n
            x = 0
        elseif k==0
            x = 1/factorial(big(n))
        else
            x = sum([X("β",(n-1-s,k-1),St)/factorial(big(s+1)) for s=0:n-1-k])
        end
        px("Computed β(n,k)=β(",n,",",k,") = ",x)
        x
    end

    # 4) Gather computation functions, same order than kinds list
    kcf = [("z",compute_z) , ("β",compute_β)]
    
    # 5) Create stock
    St = Stock(kcf,p)
    St
end

function get_z(oc,St)
    coords = oc.n
    X("z",coords,St)
end

function get_β(oc,St)
    coords = oc.n, oc.k
    X("β",coords,St)
end

function test_plot()
    St = create_stock()
    St.p["y"] = 0.52
    # y = 0.59 # transition

    # Computes
    n_max = 100
    ns = range(1, n_max)
    zs = [X("z",n,St) for n in ns]

    # Builds a plot
    f = Figure()
    ax = Axis(f[1, 1])
    scatter!(ax, ns, zs)
    return f
end

test_plot()
