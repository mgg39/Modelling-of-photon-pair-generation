function sinc = sinc(x, Zmax)
    sinc = sin(x.*Zmax/2).*(x.*Zmax/2).^-1;
end