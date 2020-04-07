function [f, Df, D2f] = fpascal()
    A = pascal(4);
    b = -ones(4,1);
    
    f   = @(x)  0.5*dot(x, A*x) + dot(b,x) +1;
    Df  = @(x)  A*x +b;
    D2f = @(x)  A;
end

