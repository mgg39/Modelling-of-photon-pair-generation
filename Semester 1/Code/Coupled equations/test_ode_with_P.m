clear all
clc

gamma=1;
kappa=0.1;
g = 1;

xini=[1,0,0];

options = odeset('RelTol',1e-8,'AbsTol',1e-10);

zmax=50;

[z,x]=ode45(@(z,x) rhs(z,x,gamma,kappa,g), [0 zmax], xini,options);

figure
plot(z,abs(x).^2)
legend('F','S','P')

%%

function  y= rhs(z,x,gamma,kappa,g)
 
 y=zeros(size(x));

 y(1)=1i*gamma*conj(x(1)).*x(2).*exp(1i*kappa*z);
 y(2)=0.5i*gamma*x(1).^2.*exp(-1i*kappa*z) - 1i*g*x(3);
 y(3)=-1i*g*x(2);
 
end