function [x, iter] = lsBFGS(f, x0, itmax)
% Purpose: approximate a local min of f using the linesearch algorithm
% and the (iBGFS) update formula (to avoid the solution of linear systems)
%
% Same parameters and results as lineDFP,
% with the exception that the (iBGFS) update formula is used.
%
    iter = 0;
    n = length(x0);
    I = speye(n);
    H = speye(n);
    gk = apGrad(f, x0);
    gnew = 1;
    x = x0;
    tol=1e-5;
    while norm(gnew) > tol
        dk = -H*gk;
        
        [alpha, gnew] = lineSearch(f, x, dk, gk);
        
        s = alpha*dk;
        x = x + s;
        gamma = gnew - gk;
        rhoinv = dot(s, gamma);
        
        H = (I - s*gamma'/rhoinv)*H*(I - gamma*s'/rhoinv) + (s*s')/rhoinv;  % O(n^3)
        gk = gnew;
        iter = iter+1;   % STOP si es grade.
        if  iter >= itmax
            break
        end
    end
end