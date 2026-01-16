# to find the two states which are at the Dirac point at K
# need to have V=0 !!!
function find_zero_energy_modes(Es;tol=1e-3,name="")
    # px("Energies ",Es)
    r = -1
    for j=1:length(Es)-1
        if sum(abs.(Es[j:j+1])) < tol # because V=0
            r = j
            break
        end
    end
    if r == -1 # if not found
        M = -1 # will be the first positive eigenvalue
        for j=1:length(Es)
            if Es[j] > 0 M = j end
        end
        # px("Energies are ",Es)
        px("ZERO ENERGY NOT FOUND FOR ",name,", PROBLEM, take 1. Should increase precision ? Around 0, energies are")
        plot_Es_around(Es,M; amplitude=4)
        r = 1
        # @assert false
    end
    # px("Find zero energy eigenmodes, r=",r," among ",length(Es)," eigenvalues")#," Es=",Es)
    # px("Search for zero energy eigenmodes :",j)
    r, Es[r]
end
#
# check that the list contains at least two positive values
function check_at_least_two_positive(Es)
    n = length(filter(x-> x>0, Es))
    if n ≤ 1
        px("PROBLEM ENERGIES : NEED TWO POSITIVE AT LEAST")
        px(Es)
        @assert false
    end
end

function plot_Es_around(Es,i_center; amplitude=-1)
    ampli = amplitude==-1 ? floor(Int,length(Es)/2) : amplitude
    m = max(1,i_center-ampli)
    M = min(length(Es),i_center+ampli)
    px("Exact energies ",Es[m:M])
end

function i_modes_stored(e,V,invε)
    d = Dict()
    # Exact
    d[("ex",0.5,9)] = 89
    d[("ex",0.5,10)] = 110
    d[("ex",0.5,11)] = 137
    d[("ex",0.5,13)] = 190
    d[("ex",0.5,15)] = 254
    # Effective
    d[("ef",0.5,9)] = 54
    d[("ef",0.5,10)] = 57
    d[("ef",0.5,11)] = 58
    d[("ef",0.5,12)] = 61
    d[("ef",0.5,15)] = 60
    req = (e,V,invε)
    stored = haskey(d,req)
    if stored
        return d(e,V,invε)
    end
    -1
end

function oc_V_becomes_0(oc)
    oc_V_is_0 = oc; oc_V_is_0.V_intensity = 0
    oc_V_is_0
end

# depends on everything but actually won't necessitates much computations for redundant data, like different kred
function compute_i_dmode(oc_list,St_sol::Stock) 
    # necessarily for V=0 because of find_zero_energy_modes
    # depends on everything but kred
    oc = list_to_oc(oc_list)
    req = (oc.e,oc.V_intensity,oc.invε)
    i_mode = i_modes_stored(req...)
    stored = i_mode != -1
    if stored
        px("i_modes was stored for request ",req)
        return i_mode
    end
    is_0 = abs(oc.V_intensity)+norm(oc.kred)<1e-14
    oc0 = oc
    oc0.V_intensity = 0
    oc0.kred = [0,0]
    if !is_0 # if V!=0, return the value of when V=0 by default
        return compute_i_dmode(oc_to_list(oc0),St_sol)
    end
    # now here V=0 and k=0
    @assert is_0
    Es = get_Es_unshifted(oc0,St_sol)
    i_mode, _ = find_zero_energy_modes(Es;name=oc.e)
    err = abs(Es[i_mode])
    if err > 1e-5
        px("ERROR !!!! on dmode energy is ",err)
        plot_Es_around(Es,i_mode; amplitude=5)
        # @assert false
    end
    i_mode
end

function get_i_dmode(oc,St_sol)
    X("i_dmode",oc_to_list(oc),St_sol)
end

