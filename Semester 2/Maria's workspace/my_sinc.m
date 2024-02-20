function y = my_sinc(x,Zmax)

    % Calculate sinc function values
    y = sin(x.*Zmax/2).*(x.*Zmax/2).^-1;
    
end
