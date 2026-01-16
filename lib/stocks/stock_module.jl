# Module which enables to manage recursive quantities computation
px = println

mutable struct Stock
    Xs # quantities to be recursively computed
    kinds # there are length(X) different recursive quantities
    computation_functions
    p # parameters
    debug
    # if together, then the elements of kinds and computation_functions come in pair in kinds_and_comp
    function Stock(kinds_and_computation_functions,parameters;debug=false,together=true,kinds="",computation_functions="")
        St = new()
        St.debug = debug
        # kinds
        kinds_, comp_f = kinds, computation_functions
        if together
            kinds_ = [v[1] for v in kinds_and_computation_functions]
            comp_f = [v[2] for v in kinds_and_computation_functions]
        end
        # First stuff
        St.kinds = kinds_
        St.p = parameters
        # X
        St.Xs = Dict()
        for k in St.kinds
            St.Xs[k] = Dict()
        end
        # Computation functions
        St.computation_functions = Dict()
        for i in 1:length(St.kinds)
            k = St.kinds[i]
            St.computation_functions[k] = comp_f[i]
        end
        St
    end
end

# Needs a value, computes it if it was not computed previously
function X(kind,coords,St::Stock)
    if !haskey(St.Xs[kind],coords)
        if St.debug
            px("Computes ",kind," ",coords)
        end
        x = St.computation_functions[kind](coords,St)
        St.Xs[kind][coords] = x
    end
    St.Xs[kind][coords]
end

function sets_value(x,kind,coords,St::Stock)
    if haskey(St.Xs[kind],coords)
        print("FOR TYPE ",kind," CANNOT SET VALUE ",x," WITH COORDS ",coords," because it is already set")
    end
    St.Xs[kind][coords] = x
end

# gives all the elements of kind kind
all_elements(kind,St::Stock) = [v for v in S.Xs[kind]]
