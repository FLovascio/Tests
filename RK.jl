module RK
export integrateRK, ButcherTab, ivpSystem, rk4

struct ButcherTab
	c::Array{Float64,1}
	b::Array{Float64,1}
	a::Array{Float64,2}
end

rk4=ButcherTab( [0.0,0.5,0.5,1.0],[0.16666666,0.333333,0.3333333,0.1666666],[0 0 0 0;0.5 0 0 0;0 0.5 0 0;0 0 1 0])

struct ivpSystem
	yDash::Array{Function}
	init::Array{Float64}
end

function stepRK(problem::Array{Function},init_t::Float64,init::Array{Float64},τ,tableau::ButcherTab,k_array)
	for (j,f) ∈ enumerate(problem)
		k_array[j,1]=τ*f(init_t, init)
	end
	for i=2:length(tableau.c)
		for (j,f) ∈ enumerate(problem)
			k_array[j,i]=τ*f(init_t+τ*tableau.c[i], init.+sum(transpose(tableau.a[i,1:(i-1)]).*k_array[:,1:(i-1)],dims=2))
		end
	end
	return init.+sum(transpose(tableau.b).*k_array,dims=2)
end

function integrateRK(problem::Array{Function},init::Array{Float64},t_range,tableau::ButcherTab)
	T=collect(t_range)
	τ=T[2]-T[1]
	solution=zeros(length(problem),length(T))
	k_array=zeros(length(problem),length(tableau.c))
	solution[:,1]=init
	for i=2:length(T)
		solution[:,i]=stepRK(problem,T[i-1],solution[:,i-1],τ,tableau,k_array)
	end
	return Dict([("t",T),("y",solution)])
end

function integrateRK(system::ivpSystem,t_range,tableau::ButcherTab)
	T=collect(t_range)
	τ=T[2]-T[1]
	problem=system.yDash
	solution=zeros(length(problem),length(T))
	k_array=zeros(length(problem),length(tableau.c))
	solution[:,1]=system.init
	for i=2:length(T)
		solution[:,i]=stepRK(problem,T[i-1],solution[:,i-1],τ,tableau,k_array)
	end
	return Dict([("t",T),("y",solution)])
end

end

