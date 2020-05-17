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
    B = eye(n);
    H = B;
    
    while norm(g) > tol && iter < itmax
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
        redm = -(dot(g, s) + 0.5 * dot(s, B * s);
        gamma = apGrad(xk + s) - g;
        
        %PASO 3
        if redf / redm > eta
            xk = xk + s;
            g = apGrad(f, xk);
            iter = iter + 1;
        end
        
        %PASO 4
        if redf / redm > 0.75
            
    
    

    end 

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
