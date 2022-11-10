using OrdinaryDiffEq
using LinearAlgebra
using Plots
using LaTeXStrings

Vd = 12;
rC = 0.12;
rL = 0.176;
R  = 50;
L  = 100e-6;
C  = 30e-6;
d  = 0.677;
dp = 1-d;
x₀ = [0; 0];

# Equilibria 

VC = (Vd*dp*R)/(rL+((dp*rC*R)/(R+rC))+(((dp^2)*(R^2))/(R+rC)));
IL = (Vd)/(rL+((dp*rC*R)/(R+rC))+(((dp^2)*(R^2))/(R+rC)));


function dxdt(dx,x,p,t)
    global L,C,Vd,rC,rL,dp;
    dx[1] = (1/L)*((-(rL+((dp*rC*R)/(R+rC)))*x[1]) - (((dp*R)/(R+rC))*x[2]) + Vd);
    dx[2] = (1/C)*((((dp*R)/(R+rC))*x[1]) - ((1/(R+rC))*x[2]));
end 

t = (0,0.01)

prob = ODEProblem(dxdt,x₀,t)

# Simulate using ODE Solver
sol = solve(prob, Tsit5())
 
# And plot
p = plot(sol,
         layout=(2,1), 
         xlabel=["Time \$t\$ [s]" "Time \$t\$ [s]"],
         ylabel=["Current [A]" "Voltage [V]"],
         label = ["\$i_L(t)\$" "\$v_C(t)\$"], 
         title = "Response of the Boost Converter",
         color = [:red :blue])

