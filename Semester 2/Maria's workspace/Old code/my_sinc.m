function y = my_sinc(x)
 
    % Calculate sinc function values
   y = sin(x)./x;
   y(x==0) = 1;
   
end