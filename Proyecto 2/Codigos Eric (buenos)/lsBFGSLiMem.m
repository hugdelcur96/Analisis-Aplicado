function [x, iter] = lsBFGSLiMem(f, x0, itmax, m)
% Purpose: use limited memory updates to iapproximate the inverse of the 
%          Hessian. (linesearch algorithm)
%
    iter = 0;
    n = length(x0);
    % m = 3;
    SS = [];
    GG = [];
    gk = apGrad(f, x0);
    gnew = 1;
    x = x0;
    tol = 1e-5;
    while norm(gnew) > tol
        dk = -evaluate_Hg(SS, GG, gk);
        
        [alpha, gnew] = lineSearch( f, x, dk, gk );
        
        s = alpha*dk;
        x = x + s;
        gamma = gnew - gk;
        %rhoinv = dot(s, gamma);
        
        % we have to refresh the memory
        if iter <= m
            SS = [s, SS];
            GG = [gamma, GG];
        else
            SS = [s, SS(:, 1:m-1)];
            GG = [gamma, GG(:, 1:m-1)];
        end
        
        gk = gnew;
        iter = iter+1;   % STOP si es grade.
        if  iter >= itmax
            break
        end
    end
end


function r = evaluate_Hg(SS, GG, gk)
    % SS ... s(k-1), s(k-2), ... s(k-m)
    % GG ... gam(k-1), gam(k-2), ... gam(k-m)
    % gk ... gradiente en xk.
    
    [~, m] = size(GG);
    
    if m == 0
        r = gk;
    else
        alphas = zeros(m,1);
        r = gk;
        for j = 1:m
            alphas(j) = dot( SS(:,j), r )/dot( SS(:,j),GG(:,j) );
            r = r - alphas(j)*GG(:,j);
        end
        
        r = dot( SS(:,1),GG(:,1) )/dot( GG(:,1),GG(:,1) )*r;
        
        for j = m:-1:1
            beta = dot( GG(:,j), r )/dot( SS(:,j),GG(:,j) );
            r = r + (alphas(j)-beta)*SS(:,j);
        end
    end
end