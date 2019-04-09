module RK
export integrateRK, ButcherTab, rk4

struct ButcherTab
	c::Array{Float64,1}
	b::Array{Float64,1}
	a::Array{Float64,2}
end

function stepRK(derivative::Function,init::Tuple{Float64,Float64},τ,tableau::ButcherTab,k_array)
	k_array[1]=τ*derivative(init[1], init[2])
	for i=2:length(k_array)
		k_array[i]=τ*derivative(init[1]+τ*tableau.c[i], init[2]+sum(tableau.a[i,1:(i-1)].*k_array[1:(i-1)]))
	end
	return init[2]+sum(k_array.*tableau.b)
end


function integrateRK(derivative::Function,init::Tuple{Float64,Float64},τ,t_range,tableau::ButcherTab)
	T=collect(t_range)
	solution=deepcopy(T)
	k_array=deepcopy(tableau.c)
	solution[1]=init[2]
	for i=2:length(solution)
		solution[i]=stepRK(derivative,(T[i-1],solution[i-1]),τ,tableau,k_array)
	end
	return (T,solution)
end

rk4=ButcherTab( [0.0,0.5,0.5,1.0],[0.16666666,0.333333,0.3333333,0.1666666],[0 0 0 0;0.5 0 0 0;0 0.5 0 0;0 0 1 0])

end

