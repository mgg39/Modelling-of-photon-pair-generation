%Equation
function duhatdt = rhspde(t,uhat,kappa,a,b)
duhatdt = -a*(kappa.^2)'.*uhat -b*1i*kappa*uhat; %linear and diagonal
end