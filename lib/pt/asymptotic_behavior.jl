# Measure an asymptotic behavior

function get_λs(n_start,n_end;res=20)
    # λs = logrange(10^float(n_start),10^float(n_end),res)
    λs = [10^float(x) for x in range(n_start,n_end,res)]
    λs
end

# Log and linear
function get_λs_log_lin(n_start,n_end;res=20)
    hres = Int(floor(res/2))
    hres = res
    λ_log = get_λs(n_start,n_end;res=hres)
    λ_lin = Array(range(10^float(n_start), 10^float(n_end), hres))
    l = unique(vcat(λ_log,λ_lin)) # unique erases the duplicates
    sort(l)
end

# comparison ∈ ["last","middle"]
function print_results(λs,results,reference,should;title="",comparison="last")
    nn = length(λs)
    abs_res = abs.(results)
    # Prints results
    lres = log.(abs_res)/log(10)
    println("----- "*title)
    # println("Distances (10^d, check it does not saturate beyond -33) : ",lres)
    # Prints results in relative distances
    rel_results = round.(log.(abs_res/reference)/log(10) ,digits=5)
    println("Relative distances (10^d) : ",rel_results)
    if nn ≥ 2
        # Asymptotic behavior
        logs = []
        middle_i = floor(Int,nn/2) # point of comparision needs to be in the perturbative zone
        for i=1:nn
            I_compare = comparison == "last" ? i-1 : middle_i
            if i!= (comparison == "last" ? 1 : middle_i)
                λ = λs[i]
                x = (lres[i] - lres[I_compare]) / (log(λ)/log(10) - log(λs[I_compare])/log(10))
                y = round(x,digits=5)
                push!(logs,y)
            end
        end
        println("Prints logs behavior ("*comparison*") (should be ",should," on the left) : ",logs)
    end
end

# the good one to use, so that we have relative distances
# function print_results_quotiented(λs,results,should)
    # nn = length(λs)
    # logs = []
    # Prints results
    # I_compare = floor(Int,nn/2) # point of comparision needs to be in the perturbative zone
    # lres = log.(results)
    # for i=1:nn
        # if i!=I_compare
            # λ = λs[i]
            # x = (lres[i] - lres[I_compare]) / (log(λ) - log(λs[I_compare]))
            # y = round(x,digits=10)
            # push!(logs,y)
        # end
    # end
    # println("Prints distances (check it does not saturate beyond -33) : ",lres)
    # println("Prints logs (should be ",should," on the left) : ",logs,"\n")
# end

function test()
    f(x) = 5*(x*(1+0.1*rand()))^3
    λs = get_λs(-5,-1;res=20)
    print_results(λs,f.(λs),3)
end

# test()
