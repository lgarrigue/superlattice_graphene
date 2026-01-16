px = println

# + of functions
import Base.+  
+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)  

f(x) = 4x
h(x) = 5x

a = [f,h]

function number_of_arguments(f) # gives the number of arguments of a function
	first(methods(f)).nargs - 1
end

function zero_function(n)
	if n == 1
		return x -> 0
	elseif n == 2
		return (x,y) -> 0
	else
		return (x,y,z) -> 0
	end
end

function sum_functions(a)
	m = number_of_arguments(a[1])
	g = zero_function(m)
	n = length(a)
	for i=1:n
		@assert number_of_arguments(a[i]) == m
		g = g + a[i]
	end
	g
end

g = sum_functions(a)
g(1)
