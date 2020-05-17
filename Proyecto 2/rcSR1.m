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
    deltaMax = 1.25;
    r = 1e-6;
    
    n = length(x0);
    delta = deltaMax;
    iter = 0;
    
    xk = x0;
    g = apGrad(f, xk);
    B = apHess(f, xk);
    H = speye(n);
    
    while norm(g, 'inf') > tol %&& iter < itmax
        % PASO 1
        s = -H * g;
        if dot(s, g) < 0
            if norm(s) > delta
                s = delta * s / norm(s);
            end
        else
            s = pCauchy(B, g, delta);
        end
        
        %PASO 2
        redf = f(xk) - f(xk + s);
        redm = -(dot(g, s) + 0.5 * dot(s, B * s));
        gamma = apGrad(f, xk + s) - g;
        
        cociente = redf / redm;
        
        %PASO 3
        if cociente > eta
            xk = xk + s;
            g = apGrad(f, xk);
            iter = iter + 1;
            if iter == itmax
                break
            end
        end
        
        %PASO 4
        if cociente > 0.75
            if norm(s) > 0.8 * delta
                %delta = 2 * delta;
                delta = min(2*delta, deltaMax);
            end
        end
        
        %PASO 5
        if cociente < 0.1
            delta = 0.5 * delta;
        end
        
        %PASO 6
        u = gamma - B * s;
        if abs(dot(u, s)) >= r * norm(s) * norm(u)
            B = B + u * u' / dot(u, s);
            v = s - H * gamma;
            H = H + v * v' / dot(v, gamma);
        end
    end
    x = xk;
end