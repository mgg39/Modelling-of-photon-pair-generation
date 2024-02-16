function sinc = sinc(x, Zmax)
    n =size(x);
    S = zeros(n(1), n(2));
    for c=1:n(1)
        for d = 1:n(2)
            if x(c, d) == 0
                S(c, d) = 1;
            else
            S(c, d) = sin(x(c, d)*Zmax/2)/(x(c, d)*Zmax/2);
            end
        end
    end
    sinc = S;
end