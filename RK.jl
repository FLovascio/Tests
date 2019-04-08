module RK
export RK

struct ButcherTab
	c::Array{Float64,1}
	b::Array{Float64,1}
	a::Array{Float64,2}
end

function RK(derivative::Function,init::Tuple{Float64},t_range,tableau::ButcherTab)
	T=collect(t_range)
	solution=deepcopy(T)
	k_array=deepcopy(tableau.c)
	for (i,t)∈enumerate(T)
		for i=1:length(k_array)

end

