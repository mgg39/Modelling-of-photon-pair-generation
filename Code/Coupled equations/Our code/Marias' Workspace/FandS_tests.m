clear all
clc

gamma=1;
kappa=0.1;

xini=[1,0, 0]; %add t

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

zmax=50;

[z,x]=ode45(@(z,x) rhs(z,x,gamma,kappa), [0 zmax], xini,options);

figure
plot(z,abs(x).^2)

%%

function  y= rhs(z,x,t,gamma,kappa)
 
 y=zeros(size(x));

 % UF = @(z,t,F,dFdt) i*(dFdt) - r2*(dFdt - dFdt(2:end-1)) - S.*conj(F).*exp(1i*k*z);
 y(1)=-(dx(1)dt) + 1i*gamma*(dx(1)dt + 1i*dx(1)dt(2:end-1)) + 1i*gamma*conj(x(1)).*x(2).*exp(1i*kappa*z);
 % US = @(z,t,S,dSdt) i*(dSdt) - s1*(dSdt) - s2*(dSdt - dSdt(2:end-1)) - (F.^2)/2.*exp(-1i*k*z);
 y(2)=-(dx(2)dt) + 1i*gamma*(dSdt) + 1i*gamma*(dx(2)dt + 1i*dx(2)dt(2:end-1)) + 0.5i*gamma*x(1).^2.*exp(-1i*kappa*z);
 
end