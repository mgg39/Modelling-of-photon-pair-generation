%Equation
function duhatdt = rhspde(t,uhat,kappa,a)
duhatdt = -a*(kappa.^2)'.*uhat; %linear and diagonal
end