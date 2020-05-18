function [x, iter] = rcSR1(f, x0, itmax)
% In : f     ... (handle) function to be optimized
%      x0    ... (vector) initial point
%      itmax ... (natural number) upper bound for number of iterations
%
%
% Out: x    ... (vector) last approximation of a stationary point
%      iter ... (natural number) number of iterations
    
    
    tol = 1e-5;
    eta = 0.1;
    deltaMax=1.25;
    r = 1e-6;
    n = length(x0);
    delta = deltaMax;
    iter = 0;
    x = x0;
    g = apGrad(f, x);
    B = apHess(f, x);
    H = speye(n);
    
    while norm(g, 'inf') > tol
        %PASO 1 
        s = -H * g;
        if dot(s, g) < 0
            if norm(s) > delta
                s = delta * s / norm(s);
            end
        else
            s = pCauchy(B, g, delta);
        end
         
        %PASO 2 
        cociente = -(f(x) - f(x+s)) / (dot(g, s) + 0.5 * dot(s, B*s));
        gnew = apGrad(f, x+s);
        gamma = gnew - g;
        
        %PASO 3
        if cociente > eta
            x = x + s;
            g = gnew;
            iter = iter + 1;
            if iter == itmax
                break
            end
        end
        
        if cociente > 0.75 %PASO 4
            if norm(s) > 0.8 * delta
                delta = min(2 * delta, deltaMax);
            end
        elseif cociente < 0.1 %PASO 5
            delta = 0.5 * delta;
        end
        
        v = gamma -B * s;
        
        %PASO 6
        if abs(dot(v,s)) >= r * norm(s) * norm(v)
            B = B + v * v' / dot(v, s);
            u = s - H * gamma;
            H = H + u * u' / dot(u, gamma);
        end
    end
end