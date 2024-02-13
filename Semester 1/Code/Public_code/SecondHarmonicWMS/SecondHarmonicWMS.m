%Second harmonic simulation of WMS using Schilt's mathematic derivation

x = -10:0.01:10; %normalized frequency x = (v-vline)/dvline
dv = 2.2; %laser frequency deviation
dvline = 2; %HWHM of absorption profile
I0 = 1; %input power
a0 = 1; %absorpance
psi = pi; %IM-FM phase
phi = 2.*psi; %Lockin phase
Iomega = 1; %low frequency ramp power
pw = 0.5; %modulation frequency power variation


m = dv./dvline; 
X = 1 - x.^2 + m.^2;
r = sqrt(X.^2 +4.*x.^2);



s1 =I0*a0*((sqrt(2)./m).*((-x).*sqrt(r+X) + sign(x) .* sqrt(r-X))./r);
s2 =I0*a0*(-4./m.^2 + (sqrt(2)./m.^2)*((r+1-x.^2).*sqrt(r+X)+2.*abs(x).*sqrt(r-X))./r);
s3 =(-I0*a0./m.^3)*(16.*x+(sqrt(2)./r).*(x.^3-3.*x.*(r+1)).*sqrt(r+X)+(sqrt(2)./r).*sign(x).*(1-3.*x.^2-3.*r).*sqrt(r-X));

s2p = Iomega.*cos(2.*psi.*s2) - pw.*dvline.*(m./2).*(cos(psi.*s1) + cos(3.*psi.*s3));
s2q = Iomega.*sin(2.*psi.*s2) - pw.*dvline.*(m./2).*(sin(psi.*s1) + sin(3.*psi.*s3));


s2phi = -(s2p*cos(phi) + s2q*sin(phi));

%s2max = I0.*a0.*(-4./m.^2 + (2./m.^2).*((m.^2+2)./sqrt(m.^2+1)));

plot(x,s2phi)