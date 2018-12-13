
#Pkg.add("")
#Pkg.add("PGFPlots")
using Plots;
#pgfplots();
#plotly();
gr();
using LaTeXStrings
#options(jupyter.plot_mimetypes = c("text/plain", "image/png" ));

#Define preshock parameters
const global p₀=2.0;
const global u₀=1.0;
const global c²=5.0;
const global f₀=2.5; #cs^2*(1-f_d), f_d =0.5
const global ρ₀=p₀/f₀;
const global D=ρ₀*u₀;
const global B=ρ₀*u₀^2+p₀;
global ts=0.1;

#Defining Differential EQN.
∂ₓp(p)=(((p*(B-p))-D*p₀*u₀)/D)/(ts*(c²-p*((B-p)/D^2)));

#Runge-Kutta 4
K₁(p,h)=h*∂ₓp(p);
K₂(p,h)=h*∂ₓp(p+K₁(p,h)/2);
K₃(p,h)=h*∂ₓp(p+K₂(p,h)/2);
K₄(p,h)=h*∂ₓp(p+K₃(p,h));
p₊₁(p,h)=0.166666666*(K₄(p,h)+2.0K₃(p,h)+2.0K₂(p,h)+K₁(p,h)) +p;

function Integrate(x::Array{Float64},pInit)
    p=pInit;
    P₁=pInit;
    for i=1:99999
        p=p₊₁(p,h);
        P₁=vcat(P₁,p);
    end
    P₋=deepcopy(P₁)
    for i=1:length(P₁)
        P₋[i]=P₁[(length(P₁)-i+1)]
    end
    P₋=map(x->2.0pInit-x,P₋);
    Pn=vcat(P₋,P₁);
    return Pn
end

#Calculating P at shock centre
#pInit=0.25*(B-sqrt(1.0+4.0*D^2*f₀)-2.0p₀);
pᵣ=14.0;
pInit=0.5*(p₀+pᵣ);

#running the RK solver
h=0.0001;
x=map(x->h*x,collect(-100000:99999));
Pn=Integrate(x,pInit)


plot(x,Pn,xlabel="x",ylabel="P")




