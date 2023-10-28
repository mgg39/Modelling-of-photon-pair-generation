% Testing the equations based on the fundamental deconstructed mathematics

% ALL CONSTANTS ---------------------------------------------------

kappa=0.1; % Taken From Andriy's Code, probably worth TEST
emax = 50; % TFAC -> max e: TEST
xini=[2,0]; % TFAC UNSURE what this is? -> may be n equation?

%TODO: some of the variable values should be constants that we can test,
%for now I am leaving them in the section bellow but this should be
%modified as we better understand their origin

% COUPLING CONSTANT
C = placeholder; %defined in section 1.3: C = (B+ - B-)/2 

% ALL VARIABLES ---------------------------------------------------
placeholder = 1; %used for the time being

%TODO : maybe make this variables appear as functions dependent on modulus
%m = {f,s,p}

% TIME
t = [0,1];

% PROPAGATION
% Propagation constant
% UNSURE: leaving it as 1 for now -> ideal fully conserved momentum scenario
beta_one[0,1,2] = [placeholder,placeholder,placeholder];
beta_two[0,1,2] = [placeholder,NaN,NaN];  % UNSURE = nxt step?

% Walk off length
%TODO : check t0 definition
y_distance = (2*(t(0)^(2))/abs(beta_two(1));
 
y_s = t(0)/(beta_one(2) - beta_one(1));
y_p = t(0)/(beta_one(3) - beta_one(1));

% Walk off relations
s_one = y_distance/y_s; %s1
p_one = y_distance/y_p; %s2

m_one = [0,1,2] = [NaN,s_one,p_one];

% Propagation relations
%TODO : make sure beta_f = beta_f_2 in the equations -> otherwise change
r_two = -beta_f/abs(beta_f_2);
s_two = -beta_s/abs(beta_f_2);
p_two = -beta_p/abs(beta_f_2);

m_two = [0,1,2] = [r_two,s_two,p_two];

% COUPLING CONSTANT
% UNSURE: leaving it as 1 for now -> ideal fully conserved momentum scenario
g = y_distance/C; %normalised coupling constant

% ------------------------------------------------------------------------
% EQUATIONS ---------------------------------------------------
% ------------------------------------------------------------------------

% F ---------------------------------------------------
F = m_two(1)*(F'') + 1i*conj(x(1)).*x(2).*exp(1i*kappa*e); %dF/de
%TODO UNSURE : HOW DO I EVEN CALL F''!!! -> F''=d2F/dt2

% S ---------------------------------------------------
S = m_one(2)*(S') + m_two(2)*(S'') + 0.5i*F.^2.*exp(-1i*kappa*e) + g*P;  %dS/de

% P ---------------------------------------------------
P = m_one(3)*(P') + m_two(3)*(P'') + g*S; %dP/de
 


