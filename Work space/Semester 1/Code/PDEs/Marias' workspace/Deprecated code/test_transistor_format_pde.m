%https://www.mathworks.com/help/matlab/math/solve-transistor-pde.html

clear all
clc

Beta_2 = 0.83;
r_2 = -Beta_2/(2*1i);

%solution mesh
x = linspace(0,1,50);
t = linspace(0,1,50);

%solve equation
m = 0;
eqn = @(x,t,u,dudx) PDE(x,t,u,dudx,r_2);
ic = @(x) IC(x);
sol = pdepe(m,eqn,ic,@BC,x,t);

surf(x,t,u)
title('Numerical Solution (50 mesh points)')
xlabel('Distance x')
ylabel('Time t')
zlabel('Solution u(x,t)')


%Equation
function [c,f] = PDE(x,t,u,dudx,r_2)

D = r_2;
f = D*dudx;

end

%Initial conditions
function u0 = IC(x)

u0 = 0;

end

%Boundary conditions
function [pl,ql,pr,qr] = BC(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;
end
