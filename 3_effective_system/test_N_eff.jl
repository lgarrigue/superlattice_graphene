px = println
d = Dict(1=>1, 2=>4, 3=>8, 4=>19, 5=>28, 6=>38, 7=>55, 8=>70, 9=>92, 10=>112, 11=>142, 12=>164, 13=>190, 14=>226)

keys = []
vals = []
for (k,v) in d
    push!(keys,k)
    push!(vals,v)
end
p = sortperm(keys)
keys = keys[p]
vals = vals[p]

for i=1:length(keys)
    k = keys[i]
    v = vals[i]
    x = v-(k)^2
    px(k," ",x)
end
