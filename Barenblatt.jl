
function barenblatt(x::Float64, t::Float64)
    A=0.16666666666666
    if A >= ((x^2)/(6*(t^0.666666)))
        return (t^(-0.333333333))*(A-((x^2)/(6*(t^0.666666))))
    else
        return 0
    end
end;

function Bl(x::Array{Float64,1},t::Float64)
    return map(u->barenblatt(u,t),x)
end;

global X=collect(-499:500);
X=map(x->0.01*x,X);


using Plots


#gr()


BB=Bl(X,0.1)
BF=Bl(X,0.5)
plot(X,Bl(X,0.5),label="t=0.5")
plot!(X,BB)

global h=0.01;
global NGH=1;

function Constant(u::Array{Float64,1},K::Float64)
    L=length(u)
    for i=1:NGH
        u[i]=K
        u[L-i]=K
    end
end

Bound(u)=Constant(u,0.0)

∂ₓ(u::Array{Float64,1},i::Int64) = (u[i+1]-u[i-1])/(2.0*h)
∂ₓₓ(u::Array{Float64,1},i::Int64) = (u[i+1]-2.0*u[i]+u[i-1])/(h^2)

Df(A::Float64) = A
F(u::Array{Float64,1},i::Int64) = Df(u[i])*∂ₓₓ(u,i) + ∂ₓ(u,i)^2

function b(j)
    J=convert(Float64,j)
    if j < 3
        B=0.33333333333333333
    else
        B=(J^2 + J -2)/(2J*(J+1))
    end
    return B
end

function w(s)
    S=convert(Float64,s)
    W=4/(S^2 + S -2)
    return W
end

function μ(j)
    J=convert(Float64,j)
    MU=((2J-1)/J)*(b(j)/b(j-1))
    return MU
end

function ν(j)
    J=convert(Float64,j)
    NU=-((J-1)/J)*(b(j)/b(j-2))
    return NU
end

μ_2(j,s)=μ(j)*w(s)

function γ_2(j,s)
    GAM=-μ_2(j,s)*(1-b(j-1))
    return GAM
end

function calc_t_parab(D::Function,D_array::Array{Float64,1},F::Array{Float64,1})
    D_array=map(x->D(x),F)
    t=maximum(D_array)
    t=0.5*abs((h^2)/t)
    return t
end

function Dot(F::Function,u::Array{Float64,1},result::Array{Float64,1})
    for i=(1+NGH):(length(u)-NGH)
        result[i]=F(u,i)
    end
    return result
end

function step_calc(τ::Float64,t_parab::Float64)
    delta=9+16*(τ/t_parab)
    #println(delta)
    return convert(Int32,ceil(0.5*sqrt(delta)-0.5))
end

function RKL_step(u::Array{Float64,1},D_array::Array{Float64,1},F::Function,D::Function,τ::Float64,temp::Array{Float64,2})
    s=step_calc(τ,calc_t_parab(D,D_array,u))
    Y0=u
    temp[1,:]=Y0
    temp[2,:]=Y0+μ_2(1,s)*τ*Dot(F,Y0,temp[4,:])
    Bound(temp[1,:])
    Bound(temp[2,:])
    for j=2:s
        i=((j)%3)+1
        temp[((j)%3)+1,:]=μ(j)*temp[((j-1)%3)+1,:] +
                          ν(j)*temp[((j-2)%3)+1,:] +
                          (1-μ(j)-ν(j))*Y0 +
                          μ_2(j,s)*τ*Dot(F,temp[((j-1)%3)+1,:],temp[5,:]) +
                          γ_2(j,s)*τ*Dot(F,Y0,temp[4,:])
        Bound(temp[i,:])
    end
    return temp[((s)%3)+1,:]
end
    
    
#function RK4_step(u::Array{Float64,1},D_array::Array{Float64,1},F::Function,D::Function,τ::Float64,temp::Array{Float64,2})

function nSteps(n::Int64,init::Array{Float64,1},DARRAY::Array{Float64,1},F::Function,Df::Function,τ::Float64,TempMatrix::Array{Float64,2})
	u=deepcopy(init)
	for i=1:n
		u=RKL_step(u,DARRAY,F,Df,τ,TempMatrix)
		#Sol=hcat(Sol,u)
	end
	return u
end    

	TempMatrix=zeros(5,1000)
	DARRAY=zeros(1000)
	u=convert(Array{Float64,1},BB)
	Sol=deepcopy(u)
	umin=deepcopy(u)
	τ=0.01
	u=nSteps(40,u,DARRAY,F,Df,τ,TempMatrix)

t=collect(1:40)
plot(X,u)
plot!(X,umin)
    


