function y = my_sinc(x)
    % Ensure that the denominator is not zero
    idx = (x ~= 0);
    
    % Calculate sinc function values
    y = sin(pi * x) ./ (pi * x);
    
    % Set the value at x=0 to 1
    y(~idx) = 1;
end
